# AnimeFace 2009

The face detector for anime/manga.
This is similar to [lbpcascade_animeface](https://github.com/nagadomi/lbpcascade_animeface), but it's more accurate and supports facial landmark detection.
I developed this library in 2008~2009.

Original website: http://anime.udp.jp/ (in Japanese)

![demo](https://raw.githubusercontent.com/nagadomi/animeface-2009/master/figure/imas.png)

Figure (c) namco


## Installation

Requirements
- Ruby
- ImageMagick
- gcc, make

### Ubuntu
```
sudo apt-get install libmagickwand-dev
sudo gem install rmagick
./build.sh
```

### Arch Linux
Install `imagemagick-full` and `ruby-rmagick` from [AUR](https://aur.archlinux.org/),
and then run `./build.sh`.

## Run sample code (Ruby)

```
cd animeface-ruby
ruby sample.rb <input image>
```
View at `${input_image}_out.png`

## Extract landmarks for a folder of images (Ruby)

First install Ruby package `progress_bar`. Then run:

```
cd animeface-ruby
ruby proc_folder.rb <input image folder> <output landmark file>
```

Each of the lines in the `<output landmark file>` is a JSON string corresponding to a file in `<input image folder>`. 
View the output file for more details.

## Create new dataset with animeface-ruby

1. Prepare images first.
2. Extract face images with `animeface-ruby/face_collector.rb`
```
face_collector.rb --src <image dir> --dest <output dir> --threshold <0.0~1.0, default: 0.2> --margin <0.0~, default: 0.1>
```
3. Delete false positive images using windows explorer or something.
4. Make annotation data from the filename (filename is formatted as `${orignal_file_name_without_extension}_${x}_${y}_${width}_${height}.png`, see [example](./animeface-ruby/face2xml.rb))
