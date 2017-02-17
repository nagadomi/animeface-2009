# AnimeFace 2009

The face detector for anime/manga.
This is similar to [lbpcascade_animeface](https://github.com/nagadomi/lbpcascade_animeface), but it's more accurate and supports facial landmark detection.

![figure/figure.png](imas.png)

I developed this library in 2008~2009.
Orignal website: http://anime.udp.jp/ (in Japanese)

## Install

Requirements
- Ruby
- ImageMgaick
- gcc, make

```
sudo apt-get install libmagickwand-dev
sudo gem install rmagick
./build.sh
```

## Run sample code (Ruby)

```
cd animeface-ruby
ruby sample.rb <input image>
```
View at `${input_image}_out.png`

## Create new dataset with animeface-ruby

1. Prepare images first.
2. Extract face images with `animeface-ruby/face_collector.rb`
```
face_collector.rb --src <image dir> --dest <output dir> --threshold <0.0~1.0, defualt: 0.2> --margin <0.0~, default: 0.1>
```
3. Delete false positive image with windows explorer or something.
4. Make annotation data from the filename (filename is formatted as `${orignal_file_name_without_extension}_${x}_${y}_${width}_${height}.png`)
