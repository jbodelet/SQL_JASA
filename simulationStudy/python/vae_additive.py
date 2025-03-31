# ADDITIVE VARIATIONAL AUTO-ENCODER
# ---------------------------------
# Fits an Additive VAE to the data in the directory `dir` and saves the estimated functions and latent variables in the directory `dir/vae_fit/`.
# The decoder architecture includes a bottleneck layer of size `K` that controls the number of estimated basis functions. This number should be chosen based on the dimensionality of the observations `p` -- in the simulations, it has been optimized for best performance.

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
import os
from pathlib import Path
import scipy.stats as scs
import argparse

# Set the seeds for reproducibility
np.random.seed(42)
tf.random.set_seed(42)

# Use all but one physical cores
num_threads = os.cpu_count()
threads_used = num_threads -1 if num_threads > 1 else 1

print(f"Using {threads_used} threads.")

tf.config.threading.set_intra_op_parallelism_threads(threads_used)
tf.config.threading.set_inter_op_parallelism_threads(threads_used)


# %%
# Create a parser for the file
parser = argparse.ArgumentParser()
parser.add_argument("-dir", help="the name of the directory", type=str)
parser.add_argument("-q", help="the number of latent variables", type=int, default=1)
parser.add_argument("-K", help="the number of estimated basis functions", type=int, default=5)
args = parser.parse_args()



# %%
# ## Create the sampling layer
class Sampling(layers.Layer):
    """Uses (z_mean, z_log_var) to sample z, the vector encoding a digit."""

    def call(self, inputs):
        z_mean, z_log_var = inputs
        batch = tf.shape(z_mean)[0]
        dim = tf.shape(z_mean)[1]
        epsilon = tf.keras.backend.random_normal(shape=(batch, dim))
        return z_mean + tf.exp(0.5 * z_log_var) * epsilon


# ## Create the variational autoencoder
class VAE(keras.Model):
    def __init__(self, encoder, decoder, K, estimate_variance=False, known_variance=1., **kwargs):
        super().__init__(**kwargs)
        self.encoder = encoder
        self.decoder = decoder
        self.K = K
        self.estimate_variance = estimate_variance
        self.known_variance = tf.Variable(initial_value=known_variance, trainable=False, dtype=tf.float32)
        
        if estimate_variance:
            # We parameterize variance via log_sigma2 to ensure positivity
            self.log_variance = tf.Variable(initial_value=tf.math.log(known_variance), trainable=True, dtype=tf.float32)

        self.total_loss_tracker = keras.metrics.Mean(name="total_loss")
        self.reconstruction_loss_tracker = keras.metrics.Mean(name="reconstruction_loss")
        self.kl_loss_tracker = keras.metrics.Mean(name="kl_loss")

    @property
    def metrics(self):
        return [
            self.total_loss_tracker,
            self.reconstruction_loss_tracker,
            self.kl_loss_tracker,
        ]

    def train_step(self, data):
        with tf.GradientTape() as tape:
            z_mean, z_log_var, z = self.encoder(data)
            reconstruction = self.decoder(z)

            # Estimate variance if needed
            if self.estimate_variance:
                var = tf.exp(self.log_variance)
            else:
                var = self.known_variance

            reconstruction_loss = 0.5 * tf.reduce_mean(
                tf.reduce_sum(
                    (data - reconstruction)**2 / var + tf.math.log(var) + tf.math.log(2. * np.pi),
                    axis=1
                )
            )

            kl_loss = -0.5 * (1 + z_log_var - tf.square(z_mean) - tf.exp(z_log_var))
            kl_loss = tf.reduce_mean(tf.reduce_sum(kl_loss, axis=1))
            total_loss = reconstruction_loss + kl_loss

        grads = tape.gradient(total_loss, self.trainable_weights)
        self.optimizer.apply_gradients(zip(grads, self.trainable_weights))

        self.total_loss_tracker.update_state(total_loss)
        self.reconstruction_loss_tracker.update_state(reconstruction_loss)
        self.kl_loss_tracker.update_state(kl_loss)
        return {
            "total_loss": self.total_loss_tracker.result(),
            "reconstruction_loss": self.reconstruction_loss_tracker.result(),
            "kl_loss": self.kl_loss_tracker.result(),
            "log_variance": self.log_variance if self.estimate_variance else self.known_variance,
        }


def generate_vae(input_shape, latent_dim, K=5, known_variance=1., estimate_variance=False):
    # ## Create the encoder network
    """ Create the encoder network; K is the dimension of the last layer. A small number helps not to overfit. This is determined by n, see the paper for references."""

    encoder_inputs = keras.Input(shape=input_shape)

    x = layers.Flatten()(encoder_inputs)
    x = tf.keras.layers.Dropout(.3)(x)
    x = layers.Dense(2*input_shape[0], activation="tanh")(x)

    z_mean = layers.Dense(100, activation="tanh")(x)
    z_mean = layers.Dense(50, activation="tanh")(x)
    z_mean = layers.Dense(latent_dim, name="z_mean", activation="linear")(z_mean)

    z_log_var = layers.Dense(100, activation="tanh")(x)
    z_log_var = layers.Dense(50, activation="tanh")(z_log_var)
    z_log_var = layers.Dense(latent_dim, name="z_log_var", activation="linear")(z_log_var)
    z = Sampling()([z_mean, z_log_var])

    encoder = keras.Model(encoder_inputs, [z_mean, z_log_var, z], name="encoder")

    # ## Create the Decoder network
    # For each of the latent variables, we create a NN. The last layer takes their sum.
    x_output_list = []
    latent_inputs = keras.Input(shape=(latent_dim,))

    # now split
    latent_split = layers.Lambda(lambda x: tf.split(x, latent_dim, axis=1))(latent_inputs)
    for i in range(latent_dim):
        x = layers.Dense(100, activation="tanh")(latent_split[i])

        x_output = layers.Dense(100, activation="tanh")(x)
        x_output = layers.Dense(50, activation="tanh")(x_output)
        x_output = layers.Dense(K, activation="tanh")(x_output)
        x_output = layers.Dense(input_shape[0], activation="linear")(x_output)
        x_output = layers.Reshape(input_shape)(x_output)
        x_output_list.append(x_output)

    x_output = layers.Add(name="list_of_g")(x_output_list)
    decoder = keras.Model(latent_inputs, x_output, name="decoder")

    vae = VAE(encoder, decoder, K=K, known_variance=known_variance, estimate_variance=estimate_variance)

    return(vae)


# %%
def data_preprocessing(path):
    # we pre-train the vae by randomizing the column orders
    data = pd.read_csv(path)
    data = np.array(data, dtype='float32')
    data = np.expand_dims(data, axis=-1)
    return data


# %%
def vae_fit(data, K=5, latent_dim=5, known_variance=1., estimate_variance=True):
    callback = tf.keras.callbacks.EarlyStopping(monitor='total_loss', patience=50)
    
    vae = generate_vae(input_shape=data.shape[1:], latent_dim=latent_dim, K=K, known_variance=known_variance, estimate_variance=estimate_variance)
    
    vae.compile(optimizer=keras.optimizers.Adam(learning_rate=0.001))
    vae.fit(data, epochs=1000, callbacks=[callback])

    # We now create a model to give as output the `latent_dim` functions evaluated for a given $z$.

    intermediate_model = keras.Model(inputs = vae.decoder.input,
                                        outputs = vae.decoder.get_layer("list_of_g").input)

    return(vae, intermediate_model)


def save_output(data, vae, intermediate_model, file):

    # save g
    dist = scs.norm(0,1)
    prob = np.arange(1,201)/201.
    z = dist.ppf(prob)
    z = np.repeat(z[:, np.newaxis], args.q, axis=1)

    y_list = intermediate_model.predict(z)

    if args.q == 1:
        y_list = [y_list]  # make it a list so the rest of the code works the same

    for k in range(len(y_list)):
        y = y_list[k].squeeze()
        zy = np.concatenate((z[:,0, np.newaxis], y), axis=1)
        df = pd.DataFrame(
            data=zy, 
            columns=["z"] + ["g"+str(i+1) for i in range(y.shape[1])])

        # create folder
        if args.q == 1:
            outdir = PATH_TO_FIT/f"g"
        else:
            outdir = PATH_TO_FIT/f"g{(k+1)}"

        outdir.mkdir(parents=True, exist_ok=True)

        # save
        df.to_csv(outdir/file, index=False)

    # save z

    z, _, _ = vae.encoder(data)
    df = pd.DataFrame(
        data = z,
        columns = ["z" + str(i+1) for i in range(z.shape[1])]
    )
    # create folder
    outdir = PATH_TO_FIT/"z"
    outdir.mkdir(parents=True, exist_ok=True)

    # save
    df.to_csv(outdir/file, index=False)
    

# %% List all files
PATH_TO_DATA = Path(f"{args.dir}/x/")
PATH_TO_FIT  = Path(f"{args.dir}/vae_fit/")

files = []
for (dirpath, dirnames, filenames) in os.walk(PATH_TO_DATA):
    files.extend(filenames)
    break

print(f"\n{len(files)} files found in  directory {args.dir}.")

for file in files:
    data = data_preprocessing(PATH_TO_DATA/file)
    # When estimate_variance is True (defaults), the model estimates the variance of the residuals, with initial values given by known_variance. Residuals are assumed to be iid Gaussian and homoscedastic.
    vae, intermediate_model = vae_fit(data,  K=args.K, latent_dim=args.q,  known_variance=1., estimate_variance=True)
    save_output(data, vae, intermediate_model, file)