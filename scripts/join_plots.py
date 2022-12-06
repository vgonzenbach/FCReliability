from PIL import Image, ImageFont, ImageDraw
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
import os, tempfile

def set_wd():
    WORKDIR = os.path.join(os.path.dirname(__file__), '..')
    os.chdir(WORKDIR)

def crop(path):
    plot = Image.open(path)
    l = min(plot.size)
    plot = plot.crop(box=(plot.width/2 - l/2, 0, plot.width/2 + l/2, l))
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix='.png')
    plot.save(tmp.name)

    return tmp.name, l

def add_border(path, l):
    img = plt.imread(path)
    fig, ax = plt.subplots(1)
    ax.axis('off')
    ax.imshow(img)
    arc = Arc((l/2, l/2), width=l-100, height=l-100, angle=0, theta1=0.0, theta2=95.5, color='salmon')
    ax.add_patch(arc)
    arc = Arc((l/2, l/2), width=l-100, height=l-100, angle=97.5, theta1=0.0, theta2=243, color='lightskyblue')
    ax.add_patch(arc)
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix='.png')
    plt.savefig(tmp.name, dpi=600, pad_inches=0, bbox_inches='tight', format='png')
    return tmp.name

def hjoin_plots(paths, outfile):

    plots = [Image.open(p) for p in paths]

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
def main():
    # Set project root here
    set_wd()
    dense_out = 'figures/dense.png'
    dense_paths = ['figures/dense_ICC_circle_manual-1.png', 
                   'figures/dense_discrim_circle_manual-1.png', 
                   'figures/dense_identify_circle_manual-1.png']
    dtemps = []
    for p in dense_paths:
        tmp, _ = crop(p)
        #tmp = add_border(tmp, l)
        dtemps.append(tmp)

    hjoin_plots(dtemps, dense_out)

    # Sparse
    sparse_out = 'figures/sparse.png'
    sparse_paths = ['figures/sparse_ICC_circle_manual-1.png', 
                    'figures/sparse_discrim_circle_manual-1.png', 
                    'figures/sparse_identify_circle_manual-1.png']
    
    stemps = []
    for p in sparse_paths:
        tmp, _ = crop(p)
        #tmp = add_border(tmp, l)
        stemps.append(tmp)

    hjoin_plots(stemps, sparse_out)

    # Together
    vjoin_plots([dense_out, sparse_out], 'figures/dense+sparse.png')
    return dtemps, stemps

if __name__ == "__main__":
    dtemps, stemps = main()