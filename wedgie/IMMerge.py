from PIL import Image
from IPython import embed
import glob
def merge_images(file1, file2):
    """Merge two images into one, displayed side by side
    :param file1: path to first image file
    :param file2: path to second image file
    :return: the merged Image object
    """
    image1 = Image.open(file1)
    image2 = Image.open(file2)

    (width1, height1) = image1.size
    (width2, height2) = image2.size

    result_height = height1 + height2
    result_width = max(width1, width2)

    result = Image.new('RGB', (result_width, result_height))
    result.paste(im=image1, box=(0, 0))
    result.paste(im=image2, box=(0, height1))
    return result

images1 = glob.glob('/Users/admin/HERASummer2017/HERATempTest/plots/2457747_timeavg_bug/*png')
images2 = glob.glob('/Users/admin/HERASummer2017/HERATempTest/plots/2457747_timeavg_nobug/*png')

for i in range(len(images1)):
    print images1[i]
    print images2[i]
    merged = merge_images(images1[i], images2[i])
    merged.save("{}.png".format(i))