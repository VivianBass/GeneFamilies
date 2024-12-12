
library(latex2exp)
library(ggplot2)

# Create some data
x <- seq(-2*pi, 2*pi, length.out=100)
y <- sin(x)

# Create a plot with LaTeX expressions
ggplot(data.frame(x=x, y=y), aes(x=x, y=y)) +
  geom_line() +
  ggtitle(TeX("$f(x) = \\sin(x)$")) +
  xlab(TeX("$x$ (radians)")) +
  ylab(TeX("$f(x)$")) +
  theme_minimal()



# Install dependencies first
install.packages("units")
install.packages("sf")
install.packages("gifski")
install.packages("transformr")

# Finally install gganimate
install.packages("gganimate")


library(ggplot2)
library(gganimate)
library(gifski) 


#------------------------------------------------------------------------------------------------
# Create data with multiple states
df <- data.frame(
  x = rep(1:5, 3),
  y = rnorm(15),
  state = rep(c("A", "B", "C"), each = 5)
)

# Create animated plot with state transitions
p <- ggplot(df, aes(x = x, y = y)) +
  geom_line() +
  geom_point() +
  transition_states(
    state,
    transition_length = 2,
    state_length = 1
  ) +
  labs(title = "State: {closest_state}")

animate(p, fps = 10)

#------------------------------------------------------------------------------------------------

# Example with customization
p <- ggplot(df, aes(x = time, y = value)) +
  geom_line() +
  geom_point() +
  transition_reveal(time) +
  # Add customizations
  labs(
    title = "Value over Time: Frame {frame} of {nframes}",
    x = "Time",
    y = "Value"
  ) +
  theme_minimal() +
  # Add shadow effect
  shadow_wake(wake_length = 0.1)

# Customize animation rendering
animate(p,
        fps = 30,          # Frames per second
        duration = 10,     # Total duration in seconds
        width = 800,       # Width in pixels
        height = 600,      # Height in pixels
        renderer = gifski_renderer()) # Use gifski for better qual