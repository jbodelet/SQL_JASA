# Define the Generator network

get_generator_std <- function(p, std){
  generator <- nn_module(
    "Generator",
    initialize = function(p, std) {
      self$fc1 <- nn_linear(1, 8)
      self$fc2 <- nn_linear(8, 12)
      self$fc3 <- nn_linear(12, p)
      self$activation <- nn_leaky_relu()
      self$std <- std
    },
    forward = function(input) {
      Z <- input
      latent <- Z[,1, drop = FALSE]
      noise <- Z[,2:ncol(Z)]
      x <- self$activation(self$fc1(latent))
      x <- self$activation(self$fc2(x))
      x <- self$fc3(x)
      return( x + self$std * noise )
    }
  )
  return( generator(p, std) )
}

get_discriminator <- function(p){
  nn_sequential(
    nn_linear(p, 16), nn_leaky_relu(),
    nn_linear(16, 32), nn_leaky_relu(), 
    nn_linear(32, 1)
  )  
}


# normally, std should not be known:
get_generator <- function(p, max_std = 2){
  generator <- nn_module(
    "Generator",
    initialize = function(p) {
      # Define linear layers and activation function
      self$fc1 <- nn_linear(1, 8)
      self$fc2 <- nn_linear(8, 12)
      self$fc3 <- nn_linear(12, p)
      self$activation <- nn_leaky_relu()
      # Register `std` as a learnable parameter with initial value 0.1
      self$std <- nn_parameter(torch_tensor(0, requires_grad = TRUE))
    },

    forward = function(input) {
      # Split the input tensor into latent vector and noise components
      Z <- input
      latent <- Z[, 1, drop = FALSE]  # Extract the first column as the latent vector
      noise <- Z[, 2:ncol(Z)]         # Remaining columns are considered as noise

      # Pass latent vector through the fully connected layers with activations
      x <- self$activation(self$fc1(latent))
      x <- self$activation(self$fc2(x))
      x <- self$fc3(x)

      # Combine the output with the scaled noise using the learnable `std`
      transf_std <- ( torch_tanh(self$std) + 1 ) * max_std / 2  # std is set between 0 and max_std
      return(x + transf_std * noise)
    }
  )
return( generator(p) )
}


