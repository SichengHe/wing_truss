mencoder mf://*.png -mf fps=5:type=png -ovc x264 -x264encopts bitrate=32000 -o output.avi
ffmpeg -i output.avi -acodec copy -vcodec copy output.mp4
