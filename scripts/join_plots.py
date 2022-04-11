from PIL import Image, ImageFont, ImageDraw

# Set project root here
PROJROOT = "/Users/vgonzenb/PennSIVE/FCReliability/"

# 
def hjoin_plots(paths, outfile):
    # Center image
    def crop(plot):
        l = min(plot.size)
        return plot.crop(box=(plot.width/2 - l/2, 0, plot.width/2 + l/2, l))

    plots = [crop(Image.open(p)) for p in paths]

    # Concatenate all images
    img = Image.new('RGB', (plots[0].width * len(plots), plots[0].height))
    for i, plot in enumerate(plots):
        img.paste(plot, (i*plot.width, 0))

    img.save(f"{outfile.replace('.png', '')}.png")
    print(f"Image saved on '{outfile.replace('.png', '')}.png")
    return None

def vjoin_plots(paths, outfile):
    # Center image
    
    plots = [Image.open(p) for p in paths]

    # Concatenate all images
    img = Image.new('RGB', (plots[0].width, plots[0].height * len(plots)))
    for i, plot in enumerate(plots):
        img.paste(plot, (0, i*plot.height))

    img.save(f"{outfile.replace('.png', '')}.png")
    print(f"Image saved on '{outfile.replace('.png', '')}.png")
    return None

# Dense
dense_out = '../figures/dense.png'
dense_paths = ['../figures/dense_ICC_circle_manual-1.png', 
               '../figures/dense_discrim_circle_manual-1.png', 
               '../figures/dense_identify_circle_manual-1.png']
hjoin_plots(dense_paths, dense_out)

# Sparse
sparse_out = '../figures/sparse.png'
sparse_paths = ['../figures/sparse_ICC_circle_manual-1.png', 
                '../figures/sparse_discrim_circle_manual-1.png', 
                '../figures/sparse_identify_circle_manual-1.png']
hjoin_plots(sparse_paths, sparse_out)

# Together
vjoin_plots([dense_out, sparse_out], '../figures/dense+sparse.png')