import math
import numpy
import cv2
import matplotlib.pyplot as plt
from tkinter import filedialog
import tkinter
import random

def image_selector():
    path = filedialog.askopenfilename()
    if path != "":
        name = ""
        i = -1
        while path[i] != "/":
            name = path[i] + name
            i -= 1
        print("Image fetched -> ", name, "!", sep = "")
        path = path.replace('/', '\\')
    else:
        print("Error Image not loaded!")
    return path

def reshape_to2D (lst, rows, cols):
    if rows * cols != len (lst):
        print ("Can not convert an array of ", len(lst), " elements to an array of size ", rows, "x", cols, sep = "")
        return
    newList = [[0 for i in range(cols)] for j in range(rows)]
    for i in range (rows):
        for j in range (cols):
            newList[i][j] = lst[(cols*i) + j]
    return newList

def reshape_to2D_image (lst, rows, cols):
    if rows * cols != len (lst):
        print ("Can not convert an array of ", len(lst), " elements to an array of size ", rows, "x", cols, sep = "")
        return
    newList = numpy.ndarray (shape=(rows,cols), dtype=numpy.uint8)
    for i in range (rows):
        for j in range (cols):
            newList[i][j] = lst[(cols*i) + j]
    return newList

def reshape_to1D (matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    newarr = [0 for j in range (rows * cols)]
    for i in range (rows):
        for j in range (cols):
            newarr[cols*i + j] = matrix[i][j]
    return newarr

def MAP (rows, cols, primeN):
    TT = [0 for i in range (rows*cols)]
    for i in range (rows*cols):
        TT[i] = ((primeN * i) % (rows * cols))
    return TT
    
def sineMap (S, B, y0):
    pi = math.pi
    I = [0 for i in range (len(S))]
    I[0] = y0
    for i in range (1, len(I)):
        I[i] = B * math.sin(pi * I[i-1])
    return I

def scramble_image (image, y0, B, primeNumber):
    r = len (image)
    c = len (image[0])
    print("Size of image:", str(r) + "x" + str(c))

    S = MAP(r, c, primeNumber)

    I = sineMap(S, B, y0)
    sortedI = list(I)
    sortedI.sort()

    J = [0 for i in range (len(I))]
    SB = list(J)
    scImg = numpy.ndarray (shape = (r*c,), dtype = numpy.uint8)

    img1D = reshape_to1D(image)
    
    map = {value: idx for idx, value in enumerate(I)}
    for i in range (len(J)):
        J[i] = map[sortedI[i]]
        SB[i] = S[J[i]]
        scImg[i] = img1D[SB[i]]
    
    return SB, scImg

def unscramble_image (image, y0, B, primeN):
    r = len (image)
    c = len (image[0])

    S = MAP(r, c, primeN)

    I = sineMap(S, B, y0)
    sortedI = list(I)
    sortedI.sort()

    J = [0 for i in range (len(I))]
    SB = list(J)
    map = {value: idx for idx, value in enumerate(I)}
    for i in range (len(J)):
        J[i] = map[sortedI[i]]
        SB[i] = S[J[i]]

    unScImg = numpy.ndarray (shape=(r*c,), dtype=numpy.uint8)

    img1D = reshape_to1D(image)

    map2 = {value: idx for idx, value in enumerate(SB)}
    for i in range (len(SB)):
        idx = map2[i]
        unScImg[i] = img1D[idx]

    unScImg = reshape_to2D_image(unScImg, r, c)

    return unScImg

# Rule 7
def calculateDNA (binNum, no_of_bits = -1):
    if no_of_bits > 0:
        reqBits = no_of_bits - len(binNum)
        newBits = '0' * reqBits
        binNum = newBits + binNum

    dna = ''
    for i in range (0, len(binNum), 2):
        bits = binNum[i:i+2]
        if (bits == '11'):
            dna += 'A'
        elif (bits == '00'):
            dna += 'T'
        elif (bits == '01'):
            dna += 'C'
        elif bits == '10':
            dna += 'G'
    
    return dna

def DNAencode(seq):
    dnaSeq = [0 for i in range (len(seq))]
    for i in range (len(seq)):
        bin8 = format (seq[i], 'b') [-8:]
        dnaSeq[i] = calculateDNA (bin8, 8)
    
    return dnaSeq

ADD = {}
ADD['AA'] = ADD['TG'] = ADD['CC'] = ADD['GT'] = 'A'
ADD['AT'] = ADD['TA'] = ADD['CG'] = ADD['GC'] = 'T'
ADD['AC'] = ADD['TT'] = ADD['CA'] = ADD['GG'] = 'C'
ADD['AG'] = ADD['TC'] = ADD['CT'] = ADD['GA'] = 'G'

SUB = {}
SUB['AA'] = SUB['TT'] = SUB['CC'] = SUB['GG'] = 'A'
SUB['AG'] = SUB['TA'] = SUB['CT'] = SUB['GT'] = 'T'
SUB['AC'] = SUB['TG'] = SUB['CA'] = SUB['GC'] = 'C'
SUB['AT'] = SUB['TC'] = SUB['CG'] = SUB['GA'] = 'G'

XOR = {}
XOR['AA'] = XOR['TT'] = XOR['CC'] = XOR['GG'] = 'A'
XOR['AT'] = XOR['TA'] = XOR['CG'] = XOR['GC'] = 'T'
XOR['AC'] = XOR['TG'] = XOR['CA'] = XOR['GT'] = 'C'
XOR['AG'] = XOR['TC'] = XOR['CT'] = XOR['GA'] = 'G'

def add (dna1, dna2):
    res = ''
    for i in range (len(dna1)):
        res += ADD[dna1[i] + dna2[i]]
    
    return res

def sub (dna1, dna2):
    res = ''
    for i in range (len(dna1)):
        res += SUB[dna1[i] + dna2[i]]
    
    return res

def xor (dna1, dna2):
    res = ''
    for i in range (len(dna1)):
        res += XOR[dna1[i] + dna2[i]]
    
    return res

def DNA_operations (dnaSB, dnaSI, J):
    dnaSeq = [0 for i in range (len(dnaSB))]
    for i in range (len (dnaSB)):
        rem = J[i] % 3
        if rem == 0:
            dnaSeq[i] = add (dnaSB[i], dnaSI[i])
        elif rem == 1:
            dnaSeq[i] = sub (dnaSB[i], dnaSI[i])
        else:
            dnaSeq[i] = xor (dnaSB[i], dnaSI[i])

    return dnaSeq

def calculateBin (dna):
    bin = ''
    for i in range (0, len(dna)):
        bit = dna[i]
        if (bit == 'A'):
            bin += '11'
        elif (bit == 'T'):
            bin += '00'
        elif (bit == 'C'):
            bin += '01'
        elif bit == 'G':
            bin += '10'
    
    return bin

def encryptImage (image):
    SB, ScI = scramble_image(image, key['y0'], key['B'], key['p'])
    dnaSB = DNAencode(SB)
    dnaScI = DNAencode(ScI)
    dnaSeq = DNA_operations (dnaSB, dnaScI, SB)

    encImg_1D = [0 for i in range (len(dnaSeq))]
    for i in range (len(encImg_1D)):
        binary = calculateBin (dnaSeq[i])
        encImg_1D[i] = int (binary, 2)      # converting binary to decimal
    
    encImg = reshape_to2D_image (encImg_1D, len(image), len(image[0]))
    if not saveImage(encImg, title = "Save As (Encrypted Image)"):
        show_img(encImg)

    return encImg

UNADD = {}
UNADD['AA'] = UNADD['TT'] = UNADD['CC'] = UNADD['GG'] = 'A'
UNADD['AT'] = UNADD['TC'] = UNADD['CG'] = UNADD['GA'] = 'T'
UNADD['AC'] = UNADD['TG'] = UNADD['CA'] = UNADD['GT'] = 'C'
UNADD['AG'] = UNADD['TA'] = UNADD['CT'] = UNADD['GC'] = 'G'

def unAdd (dna1, dnaRes):
    dna2 = ''
    for i in range (len(dna1)):
        dna2 += UNADD[dna1[i] + dnaRes[i]]
    
    return dna2

def DNA_operations_decrypt(dnaSB, dnaSeq, SB):
    dnaScI = [0 for i in range (len(dnaSB))]
    for i in range (len (dnaSB)):
        rem = SB[i] % 3
        if rem == 0:
            dnaScI[i] = unAdd (dnaSB[i], dnaSeq[i])
        elif rem == 1:
            dnaScI[i] = sub (dnaSB[i], dnaSeq[i])
        else:
            dnaScI[i] = xor (dnaSB[i], dnaSeq[i])
    
    return dnaScI

def DNAdecode (dnaSeq):
    seq = [0 for i in range (len(dnaSeq))]
    for i in range (len(seq)):
        binary = calculateBin (dnaSeq[i])
        seq[i] = int (binary, 2)
    
    return seq

def decryptImage (encImage):
    unScImg = unscramble_image (encImage, key['y0'], key['B'], key['p'])
    
    encImg_1D = reshape_to1D (unScImg)
    SB = [i for i in range(len(encImg_1D))]
    dnaSB = DNAencode(SB)

    dnaSeq = DNAencode (encImg_1D)
    
    dnaDecImg = DNA_operations_decrypt (dnaSB, dnaSeq, SB)
    decImg_1D = DNAdecode (dnaDecImg)
    decImg = reshape_to2D_image (decImg_1D, len(encImage), len(encImage[0]))

    if not saveImage(decImg, title = "Save As (Decrypted Image)"):
        show_img(decImg)

    return decImg

def saveImage (image, location = "", title = ""):
    root = tkinter.Tk()
    root.withdraw()
    root.update()

    if not location:
        files = [('PNG Image', '*.png')]
        if title:
            location = filedialog.asksaveasfilename (filetypes = files, title = title)
        else:
            location = filedialog.asksaveasfilename (filetypes = files)

    root.destroy()

    if location:
        if not location.endswith(".png"):
            location = f"{location.split('.')[0]}.png"
        cv2.imwrite(filename = location, img = image)

        return True
    else:
        print("No location provided!")
        return False

def plotHistogram (image, imageName = ""):
    img_1D = reshape_to1D (image)

    plt.figure()
    plt.hist(img_1D, 256, color = "#000000")    # Black colour
    plt.xlabel("Intensity")
    plt.ylabel("Pixel Count")
    
    if imageName:
        plt.title(imageName)

    files = [('PNG Image', '*.png')]
    location = filedialog.asksaveasfilename (filetypes = files, title = f"Save As (Histogram of {imageName})")

    if location:
        if not location.endswith(".png"):
            location = f"{location.split('.')[0]}.png"
        plt.savefig (location, dpi = 600)
    else:
        plt.show()

    plt.close()

def entropy (image):
    p = numpy.array ([(image == v).sum() for v in range(256)])

    p = p/p.sum()
    e = 0

    for i in range(len(p)):
        if p[i] != 0:
            e -= (p[i]*numpy.log2(p[i]))
    
    return e

def horizontalDistribution_allPixels (image):
    x_list = [0 for i in range (len(image)*(len(image[0])-1))]
    y_list = list(x_list)

    idx = 0
    for x in range (len(image)):
        for y in range (len(image[0])-1):
            x_list[idx] = image[x][y]
            y_list[idx] = image[x][y+1]
            idx += 1

    return x_list, y_list

def verticalDistribution_allPixels (image):
    x_list = [0 for i in range ((len(image)-1)*(len(image[0])))]
    y_list = list(x_list)

    idx = 0
    for x in range (len(image[0])-1):
        for y in range (len(image[0])):
            x_list[idx] = image[x][y]
            y_list[idx] = image[x+1][y]
            idx += 1

    return x_list, y_list

def diagonalDistribution_allPixels (image):
    x_list = [0 for i in range ((len(image)-1)*(len(image[0])-1))]
    y_list = list(x_list)
    
    idx = 0
    for x in range (len(image[0])-1):
        for y in range (len(image[0])-1):
            x_list[idx] = image[x][y]
            y_list[idx] = image[x+1][y+1]
            idx += 1

    return x_list, y_list

def horizontalCorrelation (image):
    x_list, y_list = horizontalDistribution_allPixels(image)
    meanX = 0
    meanY = 0
    for i in range (len(x_list)):
        meanX += x_list[i]
        meanY += y_list[i]
    meanX /= len(x_list)
    meanY /= len(y_list)

    numSum = 0
    denS1 = 0
    denS2 = 0
    for i in range (len(x_list)):
        numSum += (x_list[i] - meanX) * (y_list[i] - meanY)
        denS1 += (x_list[i] - meanX)**2
        denS2 += (y_list[i] - meanY)**2

    den = math.sqrt (denS1 * denS2)

    return numSum / den

def verticalCorrelation (image):
    x_list, y_list = verticalDistribution_allPixels(image)
    meanX = 0
    meanY = 0
    for i in range (len(x_list)):
        meanX += x_list[i]
        meanY += y_list[i]
    meanX /= len(x_list)
    meanY /= len(y_list)

    numSum = 0
    denS1 = 0
    denS2 = 0
    for i in range (len(x_list)):
        numSum += (x_list[i] - meanX) * (y_list[i] - meanY)
        denS1 += (x_list[i] - meanX)**2
        denS2 += (y_list[i] - meanY)**2

    den = math.sqrt (denS1 * denS2)

    return numSum / den

def diagonalCorrelation (image):
    x_list, y_list = diagonalDistribution_allPixels(image)
    meanX = 0
    meanY = 0
    for i in range (len(x_list)):
        meanX += x_list[i]
        meanY += y_list[i]
    meanX /= len(x_list)
    meanY /= len(y_list)

    numSum = 0
    denS1 = 0
    denS2 = 0
    for i in range (len(x_list)):
        numSum += (x_list[i] - meanX) * (y_list[i] - meanY)
        denS1 += (x_list[i] - meanX)**2
        denS2 += (y_list[i] - meanY)**2

    den = math.sqrt (denS1 * denS2)

    return numSum / den

size = 40000
x = [0 for i in range (size)]
y = list(x)

for i in range (size):
    x[i] = random.randrange (0,255)
    y[i] = random.randrange (0,255)

def horizontalDistribution (image, imageName = ""):
    x_list = [0 for i in range (size)]
    y_list = list(x_list)

    for idx in range (size):
        x_list[idx] = image[x[idx]][y[idx]]
        y_list[idx] = image[x[idx]][y[idx]+1]
    
    plt.figure()
    plt.plot (x_list, y_list, 'bo', markersize = 1)

    if imageName:
        plt.title(imageName + " - Horizontal")

    files = [('PNG Image', '*.png')]
    location = filedialog.asksaveasfilename (filetypes = files, title = f"Save As (Horizontal Distribution of {imageName})")
    if location:
        if not location.endswith(".png"):
            location = f"{location.split('.')[0]}.png"
        plt.savefig (location, dpi = 600)
    else:
        plt.show()
    plt.close()
    
def verticalDistribution (image, imageName = ""):
    x_list = [0 for i in range (size)]
    y_list = list(x_list)

    for idx in range (size):
        x_list[idx] = image[x[idx]][y[idx]]
        y_list[idx] = image[x[idx]+1][y[idx]]

    plt.figure()
    plt.plot (x_list, y_list, 'bo', markersize = 1)

    if imageName:
        plt.title(imageName + " - Vertical")

    files = [('PNG Image', '*.png')]
    location = filedialog.asksaveasfilename (filetypes = files, title = f"Save As (Vertical Distribution of {imageName})")
    if location:
        if not location.endswith(".png"):
            location = f"{location.split('.')[0]}.png"
        plt.savefig (location, dpi = 600)
    else:
        plt.show()
    plt.close()

def diagonalDistribution (image, imageName = ""):
    x_list = [0 for i in range (size)]
    y_list = list(x_list)
    
    for idx in range (size):
        x_list[idx] = image[x[idx]][y[idx]]
        y_list[idx] = image[x[idx]+1][y[idx]+1]

    plt.figure()
    plt.plot (x_list, y_list, 'bo', markersize = 1)

    if imageName:
        plt.title(imageName + " - Diagonal")

    files = [('PNG Image', '*.png')]
    location = filedialog.asksaveasfilename (filetypes = files, title = f"Save As (Diagonal Distribution of {imageName})")
    if location:
        if not location.endswith(".png"):
            location = f"{location.split('.')[0]}.png"
        plt.savefig (location, dpi = 600)
    else:
        plt.show()
    plt.close()

def show_img (image, title = "Image"):
    cv2.imshow(title, image)
    cv2.waitKey(0)

def main():
    img = cv2.imread (image_selector(), cv2.IMREAD_GRAYSCALE)

    # Encryption and Decryption
    encImg = encryptImage (img)
    decryptImage (encImg)

    # Histograms
    plotHistogram(img, imageName = "Original image")
    plotHistogram(encImg, imageName = "Encrypted image")

    # Image Entropy
    print("Entropy of Original Image =", entropy (img))
    print("Entropy of Encrypted Image =", entropy (encImg))
    print()
    print()

    # Correlations
    print("Horizontal Correlation of Original image =", horizontalCorrelation(img))
    print("Vertical Correlation of Original image =", verticalCorrelation(img))
    print("Diagonal Correlation of Original image =", diagonalCorrelation(img))
    print()
    print("Horizontal Correlation of Encrypted image =", horizontalCorrelation(encImg))
    print("Vertical Correlation of Encrypted image =", verticalCorrelation(encImg))
    print("Diagonal Correlation of Encrypted image =", diagonalCorrelation(encImg))

    # Distributions
    horizontalDistribution (img, imageName = "Original Image")
    verticalDistribution (img, imageName = "Original Image")
    diagonalDistribution (img, imageName = "Original Image")
    horizontalDistribution (encImg, imageName = "Encrypted Image")
    verticalDistribution (encImg, imageName = "Encrypted Image")
    diagonalDistribution (encImg, imageName = "Encrypted Image")

if __name__ == '__main__':
    key = {'p' : 97, 'y0' : 0.14584390650801455, 'B' : 3.8457838656}
    main ()
