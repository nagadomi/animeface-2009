require 'rexml/document'
require "optparse"

params = ARGV.getopts("", "face:", "src:", "dest:")

if params["face"].nil? || params["src"].nil?
  warn "usage: #{$0} --face <face image dir> --src <image dir> --dest out.xml "
  exit(-1)
end
dest = params["dest"] || "out.xml"

src_db = {}
Dir.entries(params["src"]).each do |file|
  path = File.join(params["src"], file)
  if File.file?(path)
    basename = File.basename(file).split(".").first
    src_db[basename] = path
  end
end
doc = REXML::Document.new
doc << REXML::XMLDecl.new('1.0', 'UTF-8')
root = doc.add_element("dataset")
Dir.entries(params["face"]).each do |file|
  path = File.join(params["face"], file)
  next unless File.file?(path)
  basename = File.basename(file, ".png")
  name_parts = basename.split("_")
  h = name_parts.pop
  w = name_parts.pop
  y = name_parts.pop
  x = name_parts.pop
  src_name = name_parts.join("_")
  face = root.add_element("face")
  face.add_element("face").add_text(path)
  face.add_element("src").add_text(src_db[src_name])
  face.add_element("x").add_text(x.to_s)
  face.add_element("y").add_text(y.to_s)
  face.add_element("width").add_text(w.to_s)
  face.add_element("height").add_text(h.to_s)
end
File.open(dest, "w") do |f|
  doc.write(f, 2)
end
