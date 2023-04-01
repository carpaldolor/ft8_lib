#This is a fork of https://github.com/kgoba/ft8_lib by Karlis Goba (YL3JG)

Added stream_ft8 for realtime decoding of FT8 from 

Live streaming from a RSP1 SDR on the 20m band

Example

miri_fm  -f14074000 -Musb  -g6  -s48000 -w 5000000  -m336  - | stream_ft8

Live streaming from a RTL_SDR on the 10m band

Example

rtl_fm -f28074000 -s48000 -g65 -Musb -l20 -p22 - | stream_ft8

Sample output

.........................................
230401_011815    00.000 Rx FT8   +34 +0.7 2047 KQ4CSU K7KSH DM65
230401_011815    00.000 Rx FT8   +31 +0.8 1703 EA7LZ WG7W DM34
230401_011815    00.000 Rx FT8   +29 +0.7  962 W8SEE KC0BSP 73
230401_011815    00.000 Rx FT8   +28 +0.8 1438 CQ VOTA N6PE DM13
230401_011815    00.000 Rx FT8   +27 +0.6  666 KE4PMP NE5DX 73
230401_011815    00.000 Rx FT8   +26 +0.7 1250 IS0KNG NE5SD -17
Decoded 14 messages
.........................................


Fork of https://github.com/kgoba/ft8_lib by Karlis Goba (YL3JG)
