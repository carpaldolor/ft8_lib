#This is a fork of https://github.com/kgoba/ft8_lib by Karlis Goba (YL3JG)

Added stream_ft8 for realtime decoding of FT8 from 

Example streaming from a RSP1 SDR on the 20m band

miri_fm  -f14074000 -Musb  -g6  -s48000 -w 5000000  -m336  - | stream_ft8

Example streaming from a RTL_SDR on the 10m band

rtl_fm -f28074000 -s48000 -g65 -Musb -l20 -p22 - | stream_ft8



# References and credits

Fork of https://github.com/kgoba/ft8_lib by Karlis Goba (YL3JG)
