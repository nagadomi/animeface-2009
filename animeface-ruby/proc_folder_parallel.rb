require "set"
require 'json'
require "rmagick"
require_relative "AnimeFace"
require "progress_bar"
require "parallel"

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

infoss = Parallel.map(files, progress: "Processing images...") do |file|
  image = Magick::ImageList.new(file)
  faces = AnimeFace::detect(image)

  res = faces.map do |ctx|
    face = ctx["face"]
    left_eye = ctx["eyes"]["left"]
    right_eye = ctx["eyes"]["right"]
    mouth = ctx["mouth"]
    likelihood = ctx["likelihood"]

    info = {
      "file" => "" + file,   # clone string object
      "face" => get_box(face),
      "left_eye" => get_box(left_eye),
      "right_eye" => get_box(right_eye),
      "mouth" => get_box(mouth),
      "likelihood" => likelihood
    }
    info
  end
end

fout = File.open(output_file, "w")
infoss.flatten(1).each do |info|
  fout << info.to_json
  fout << "\n"
end
