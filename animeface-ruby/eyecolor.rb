# coding: utf-8
# イラストの顔の目の色をクルクルするGIFアニメを作るコマンド
# AnimeFace-Ruby.tar.gzの最新版(2011/4/11 21時以降)が必要
# 
# ruby eyecolor.rb homu.jpg homu.gif

require "rubygems"
require "rmagick"
require_relative "AnimeFace"

module AnimeFace
  module EyeColor
    def self.change(src, face, hsl, dest = nil)
      dest ||= src.copy
      change_eyecolor(dest, src, eye_info(face["eyes"]["left"]), hsl)
      change_eyecolor(dest, src, eye_info(face["eyes"]["right"]), hsl)
      dest
    end
    
    private
    S = Magick::QuantumRange * 0.5
    B = 0.6
    A = 1.0
    
    def self.change_eyecolor(dest, src, eye, hsl)
      rect = {
        :x1 => [eye[:center][:x] - eye[:r] * 1.5, 0].max.to_i,
        :y1 => [eye[:center][:y] - eye[:r] * 1.5, 0].max.to_i,
        :x2 => [eye[:center][:x] + eye[:r] * 1.5, src.columns].min.to_i,
        :y2 => [eye[:center][:y] + eye[:r] * 1.5, src.rows].min.to_i
      }
      (rect[:x1] ... rect[:x2]).each do |x|
        (rect[:y1] ... rect[:y2]).each do |y|
          p = src.pixel_color(x, y)
          # 目の中心からの距離
          point_dist = Math.sqrt((x - eye[:center][:x]) * (x - eye[:center][:x]) +
                                 (y - eye[:center][:y]) * (y - eye[:center][:y]))
          # 目の色からの距離
          color_dist = nn_color_dist(p, eye[:color][0 ... 2])
          # 重みを計算して色相を回す
          dw = gauss(point_dist, eye[:r] * B)
          cw = gauss(color_dist, S)
          dest.pixel_color(x, y, mulp(p, hsl, [(dw * cw) * A, 1.0].min))
        end
      end
    end
    def self.gauss(d, sigma)
      Math.exp(-(d * d) / (2.0 * sigma * sigma))
    end
    def self.mulp(p1, hsl2, w)
      hsl1 =  p1.to_HSL
      hue = Math.atan2((1.0 - w) * Math.sin(hsl1[0] * Math::PI * 2.0) +
                       w * Math.sin(hsl2[0] * Math::PI * 2.0),
                       (1.0 - w) * Math.cos(hsl1[0] * Math::PI * 2.0) +
                       w * Math.cos(hsl2[0] * Math::PI * 2.0)) / (Math::PI * 2.0)
      if (hue < 0)
        hue = 1.0 + hue
      end
      #saturation = ((1.0 - w) * hsl1[1] +  w * hsl2[1])
      saturation = hsl1[1]
      Magick::Pixel.from_HSL([hue, saturation, hsl1[2]])
    end
    def self.nn_color_dist(x, data)
      data.reduce(Magick::QuantumRange * 2) do |min_dist, v|
        [min_dist,
         Math.sqrt((x.red - v.red) ** 2 +
                   (x.green - v.green) ** 2 +
                   (x.blue - v.blue) ** 2)
        ].min
      end
    end
    def self.eye_info(src)
      dest = {
        :x => src["x"], :y => src["y"],
        :w => src["width"], :h => src["height"]
      }
      dest[:xr] = src["width"] / 2
      dest[:yr] = src["height"] / 2
      dest[:r] = [dest[:xr], dest[:yr]].max
      dest[:center] = {
        :x => src["x"] + dest[:xr],
        :y => src["y"] + dest[:yr]
      }
      dest[:color] = src["colors"]
      dest
    end
  end
end

def main(src, dest)
  dest ||= "anime.gif"
  image = Magick::ImageList.new(src)
  gif = Magick::ImageList.new
  
  results = AnimeFace.detect(image)
  hue = results.map{ rand }
  bg = image.first.copy
  bg.delay = 0
  gif << bg
  
  10.times do |i|
    layer = Magick::Image.new(image.columns, image.rows) do
      self.background_color = "none"
    end
    results.each_index do |j|
      hue[j] = hue[j] > 0.9 ? hue[j] + 0.1 - 1.0 : hue[j] + 0.1
      AnimeFace::EyeColor.change(bg, results[j],
                                 [hue[j], 0.3, 1.0], layer)
    end
    layer.delay = (Math.exp((-(i % 5) ** 2 ) / (2.0 * 4.0 ** 2)) * 10).to_i
    gif << layer
  end
  gif.write(dest)
end
if (ARGV[0])
  main(ARGV[0], ARGV[1])
else
  warn "usage: #{$0} src [dest]"
  exit(-1)
end
