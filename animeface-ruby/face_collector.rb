require "pp"
require "rmagick"
require_relative "AnimeFace"
require "optparse"
require "fileutils"
require "json"

def rotate_image(image, ctx)
  left_eye = ctx["eyes"]["left"]
  lex = left_eye["x"]
  lew = left_eye["width"]
  ley = left_eye["y"]
  leh = left_eye["height"]

  right_eye = ctx["eyes"]["right"]
  rex = right_eye["x"]
  rew = right_eye["width"]
  rey = right_eye["y"]
  reh = right_eye["height"]

  center_eye_x = ((lex + lew / 2) + (rex + rew / 2)) / 2
  center_eye_y = ((ley + leh / 2) + (rey + reh / 2)) / 2

  chin = ctx["chin"]
  cx = chin["x"]
  cy = chin["y"]
  cw = chin["width"]
  ch = chin["height"]

  radian = Math.atan2(cy - center_eye_y, cx - center_eye_x)
  degree = radian * 180 / Math::PI
  image.rotate(90 - degree)
end

params = ARGV.getopts("", "src:", "dest:", "threshold:", "margin:", "rotate")

if params["src"].nil? || params["dest"].nil?
  warn "usage: #{$0} --src <image dir> --dest <output dir> --threshold <0.0~1.0, default: 0.2> --margin <0.0~, default: 0.1>"
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
        y2 = [y + (face["height"] + face["height"] * margin * 2).to_i, image.rows].min
        
        if x2 - x != y2 - y
          w = [x2 - x, y2 -y].min
          x2 = x + w
          y2 = y + w
        end

        image = rotate_image(image, ctx) if params["rotate"]

        crop = image.crop(x, y, x2 - x, y2 - y, true)
        crop.write(File.join(params["dest"], 
                             sprintf("%s_%d_%d_%d_%d.png", 
                                     File.basename(file).split(".").first,
                                     x, y, x2 - x, y2 - y)))
        crop.dispose
      end
      image.dispose
    rescue => e
      warn e.message
    end
  end
end
