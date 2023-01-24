import cv2
import os
import numpy as np
import timeit


def read_by_line(filename):
    file1 = open(filename, 'r')
    lines = file1.readlines()
    return lines


def load_images_from_folder(folder):
    image_array = []
    filenames = []
    for filename in os.listdir(folder):
        img = cv2.imread(os.path.join(folder, filename))
        if img is not None:
            image_array.append(img)
            filenames.append(filename)
    return image_array, filenames


# if an array has zero add 1 to get rid of 0 divisions
def normalize_values(val):
    val[val == 0] = 1
    val /= np.sum(val)
    return val


# splits image:im into (n_of_grid x n_of_grid) grids
def split_image(im, n_of_grid):
    row = im.shape[0] // n_of_grid
    col = im.shape[1] // n_of_grid
    pieces = [im[i:(i+row), j:(j+col)] for i in range(0, im.shape[0], row) for j in range(0, im.shape[1], col)]
    return pieces


def create_histogram(image, n_bins):
    r = image[:, :, 0]
    g = image[:, :, 1]
    b = image[:, :, 2]
    red_freq = np.zeros(n_bins)
    green_freq = np.zeros(n_bins)
    blue_freq = np.zeros(n_bins)
    step_size = 256 / n_bins
    for i in range(n_bins):
        red_freq[i] = len(r[np.where(np.logical_and(r >= i * step_size, r < (i + 1) * step_size))[0]])
        green_freq[i] = len(g[np.where(np.logical_and(g >= i * step_size, g < (i + 1) * step_size))[0]])
        blue_freq[i] = len(b[np.where(np.logical_and(b >= i * step_size, b < (i + 1) * step_size))[0]])
    return red_freq, green_freq, blue_freq


def histogram_3d(image, n_bins):
    histogram = np.zeros(n_bins**3)
    red_freq, green_freq, blue_freq = create_histogram(image, n_bins)
    count = 0
    for r in red_freq:
        for g in green_freq:
            for b in blue_freq:
                histogram[count] = r + g + b
                count += 1
    return histogram


def kld(q, s):
    return np.sum(q * np.log(q / s))


def jsd(q, s):
    m = (q + s) / 2
    return (kld(q, m) + kld(s, m)) / 2


def per_channel_avg(red1, green1, blue1, red2, green2, blue2):
    red1, green1, blue1 = normalize_values(red1), normalize_values(green1), normalize_values(blue1)
    red2, green2, blue2 = normalize_values(red2), normalize_values(green2), normalize_values(blue2)
    total = jsd(red1, red2) + jsd(green1, green2) + jsd(blue1, blue2)
    return total / 3


def calculate_divergence(image1, image2, number_of_bins, per_channel=True):
    if per_channel:
        red1, green1, blue1 = create_histogram(image1, number_of_bins)
        red2, green2, blue2 = create_histogram(image2, number_of_bins)
        return per_channel_avg(red1, green1, blue1, red2, green2, blue2)
    else:
        hist1 = histogram_3d(image1, number_of_bins)
        hist2 = histogram_3d(image2, number_of_bins)
        hist1, hist2 = normalize_values(hist1), normalize_values(hist2)
        return jsd(hist1, hist2)


def calculate_divergence_with_grids(im1, im2, n_of_bins, n_grids, per_channel=True):
    split1 = split_image(im1, n_grids)
    split2 = split_image(im2, n_grids)
    total = 0
    for s1, s2 in zip(split1, split2):
        total += calculate_divergence(s1, s2, n_of_bins, per_channel)
    total /= len(split1)
    return total


def instance_recognition(query_instance, support_set, n_bins, n_grid=0, per_channel=True):
    min_index = -1
    min_value = float('inf')
    if n_grid == 0:
        for idx, support in enumerate(support_set):
            value = calculate_divergence(query_instance, support, n_bins, per_channel)
            if value < min_value:
                min_index = idx
                min_value = value
    else:
        for idx, support in enumerate(support_set):
            value = calculate_divergence_with_grids(query_instance, support, n_bins, n_grid,  per_channel)
            if value < min_value:
                min_index = idx
                min_value = value
    return min_index


def calculate_accuracy(query_set, support_set, query_filename, support_filename, n_bins, n_grid=0, per_channel=True):
    correct = 0
    for idx, query in enumerate(query_set):
        index = instance_recognition(query, support_set, n_bins, n_grid, per_channel)
        if query_filename[idx] == support_filename[index]:
            correct += 1
    return (correct / len(query_set)) * 100


if __name__ == '__main__':
    query1, filename_query1 = load_images_from_folder('dataset/query_1')
    query2, filename_query2 = load_images_from_folder('dataset/query_2')
    query3, filename_query3 = load_images_from_folder('dataset/query_3')
    images_support, filename_support = load_images_from_folder('dataset/support_96')

    # 3d n_bins = 8
    # per_channel n_bins = 16

    start = timeit.default_timer()
    accuracy = calculate_accuracy(query1, images_support, filename_query1, filename_support, 2, False)
    print('Accuracy: {}%'.format(accuracy))
    stop = timeit.default_timer()
    print('Time: ', stop - start)

    start = timeit.default_timer()
    accuracy = calculate_accuracy(query1, images_support, filename_query1, filename_support, 2, True)
    print('Accuracy: {}%'.format(accuracy))
    stop = timeit.default_timer()
    print('Time: ', stop - start)


