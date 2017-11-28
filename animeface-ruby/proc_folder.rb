require "set"
require 'json'
require "rmagick"
require_relative "AnimeFace"
require "progress_bar"

def get_box(thing)
  keys = Set["x", "y", "width", "height"]
  return thing.select { |key, value| keys.include? key }
end

if ARGV.size != 2
  warn "Usage: #{$0} <dir_path> <output_file>"
  exit(-1)
end

dir_path = ARGV[0]
files = Dir["#{dir_path}/*"]

output_file = ARGV[1]
fout = File.open(output_file, "w")

bar = ProgressBar.new(files.length)
files.each do |file|

  image = Magick::ImageList.new(file)
  faces = AnimeFace::detect(image)
  image.destroy!

  faces.each_with_index do |ctx, index|
    face = ctx["face"]
    left_eye = ctx["eyes"]["left"]
    right_eye = ctx["eyes"]["right"]
    mouth = ctx["mouth"]
    likelihood = ctx["likelihood"]

    info = {
      "file" => file,
      "face" => get_box(face),
      "left_eye" => get_box(left_eye),
      "right_eye" => get_box(right_eye),
      "mouth" => get_box(mouth),
      "likelihood" => likelihood
    }

    fout << info.to_json
    fout << "\n"
    bar.increment!
  end
end
