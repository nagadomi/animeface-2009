require "pp"
require "rmagick"
require_relative "AnimeFace"
require "optparse"
require "fileutils"
require "json"

params = ARGV.getopts("", "src:", "dest:", "threshold:")

if params["src"].nil? || params["dest"].nil?
  warn "usage: #{$0} --src <image dir> --dest <output dir> --threshold <0.0~1.0, defualt: 0.2> --margin <margin, default: 0.1>"
  exit(-1)
end
FileUtils.mkdir_p(params["dest"])
threshold = params["threshold"] ? params["threshold"].to_f : 0.2
margin = params["margin"] ? params["margin"].to_f : 0.1

Dir.entries(params["src"]).each do |file|
  if file =~ /\.(jpg|png|jpeg)$/i 
    begin
      image = Magick::ImageList.new(File.join(params["src"], file))
      faces = AnimeFace::detect(image, {:threshold => threshold})
      faces.each do |ctx|
        face = ctx["face"]
        x = ([face["x"] - face["width"] * margin, 0].max).to_i
        y = ([face["y"] - face["height"] * margin, 0].max).to_i
        x2 = [x + (face["width"] + face["width"] * margin * 2).to_i, image.columns].min
        y2 = [y + (face["width"] + face["width"] * margin * 2).to_i, image.rows].min
        
        if x2 - x != y2 - y
          w = [x2 - x, y2 -y].min
          x2 = x + w
          y2 = y + y
        end
        crop = image.crop(x, y, x2 - x, y2 - y, true)
        crop.write(File.join(params["dest"], 
                             sprintf("%s_%d_%d_%d_%d.png", 
                                     File.basename(file).split(".").first,
                                     face["x"], face["y"], face["width"], face["height"])))
      end
    rescue => e
      warn e.message
    end
  end
end
