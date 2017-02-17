require "pp"
require "rubygems"
require "RMagick"
require "AnimeFace"

image = Magick::ImageList.new("test.jpg")

result = AnimeFace::detect(image)
pp result

result = AnimeFace::detect(image, { :step => 2.0, :min_window_size => 32, :scale_factor => 1.1 })
pp result
