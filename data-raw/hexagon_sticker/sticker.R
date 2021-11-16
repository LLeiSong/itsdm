# The code to make sticker
library(sysfonts)
library(hexSticker)

imgurl <- 'data-raw/hexagon_sticker/sticker_img.png'
font_add_google("Zen Dots") # Righteous, Zen Dots, Limelight
sticker(imgurl, package="itsdm",
        p_size = 12, p_y = 1.1,
        p_family = 'Zen Dots', p_fontface = 'bold',
        s_x = 1.02, s_y = 1.05, s_width = .6,
        h_fill = 'orange', h_color = '#53222e',
        filename = "man/figures/hexagon_sticker.png",
        dpi = 600)
