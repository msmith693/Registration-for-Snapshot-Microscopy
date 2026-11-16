import numpy as np
from tqdm import tqdm
from itertools import combinations, permutations


def is_near(num, target, tolerance):
    """
    Helper function to determine if number lies within tolerance of a target number.

    :returns bool: True if number within tolerance of target.
    """
    ul = target + tolerance
    ll = target - tolerance
    return (num <= ul) and (num >= ll)


def remove_near_duplicates(data, tolerance_z, tolerance_y, tolerance_x):
    """
    Removes points within a per-axis distance specified by z, y and x values. Used to clean up localised beads, where clumping has occurred.

    :param tolerance_z: Distance in pixels for two points to be deemed distinct in the z direction.
    :param tolerance_y: Distance in pixels for two points to be deemed distinct in the y direction.
    :param tolerance_x: Distance in pixels for two points to be deemed distinct in the x direction.

    """
    results = [data[0, :]]
    for i in range(data.shape[0] - 1):
        keep = True
        for result in results:
            if (
                is_near(data[i + 1, 0], result[0], tolerance_z)
                & is_near(data[i + 1, 1], result[1], tolerance_y)
                & is_near(data[i + 1, 2], result[2], tolerance_x)
            ):
                keep = False
                break
        if keep:
            results.append(data[i + 1, :])
    return np.array(results)


def distance(coord1, coord2):
    """
    :param coord1: array of zyx coordinate points.
    :param coord2: array of zyx coordinate points.

    :return float: Euclidean distance between point A and point B in 3D space.

    """
    return np.sqrt(
        (coord1[0] - coord2[0]) ** 2
        + (coord1[1] - coord2[1]) ** 2
        + (coord1[2] - coord2[2]) ** 2
    )


def get_descriptors(bead_coordinates, threshold):
    """
    Returns all possible tetrahedral descriptors and their normalised side lengths between given coordinates under a distance threshold.

    :param bead_coordinates: Array of localised bead coordinates in zyx order.
    :param threshold: Sum of the descriptor sides. Hyperparameter to prevent generating large descriptors.

    :returns descriptors_accepted: ndarray of accepted descriptors.
    :returns descriptor_side_lengths: normalised array of corresponding accepted descriptor's side length, sorted from smallest to largest.
    """
    descriptors = list(
        combinations(bead_coordinates, 4)
    )  # for a coord set of 20, the number of triangles is 1,140 and number of descriptors is 4,845...
    descriptors_accepted = []
    descriptor_side_lengths = []
    for descriptor in tqdm(descriptors, desc="Identifying descriptors"):
        side_lengths = []
        side_lengths.append(distance(descriptor[0], descriptor[1]))
        side_lengths.append(distance(descriptor[0], descriptor[2]))
        side_lengths.append(distance(descriptor[0], descriptor[3]))
        side_lengths.append(distance(descriptor[1], descriptor[2]))
        side_lengths.append(distance(descriptor[1], descriptor[3]))
        side_lengths.append(distance(descriptor[2], descriptor[3]))
        sorted_side_lengths = np.sort(side_lengths)
        if sum(sorted_side_lengths) <= threshold:
            descriptors_accepted.append(descriptor)
            normalised_side_lengths = sorted_side_lengths / sorted_side_lengths[0]
            descriptor_side_lengths.append(normalised_side_lengths)
    descriptor_side_lengths = np.array(descriptor_side_lengths)
    return descriptors_accepted, descriptor_side_lengths


def get_abs_difference(side_lengths_1, side_lengths_2):
    """
    Calculates the sum of absolute differences in side lengths of two tetrahedral descriptors.

    :param side_lengths_1: An array containing the vertice corner points for a tetrahedron
    :param side_lengths_2: An array containing the vertice corner points for a tetrahedron
    """
    return np.sum(np.abs(side_lengths_1 - side_lengths_2))


def descriptors_close(descriptor_1, descriptor_2, tolerance):
    """
    Checks if descriptor in first volume is in a similar point in space to descriptor in second volume.
    Key function for improving algorithm's performance. Relies on fact that transform is near identity.

    :param descriptor_1: Descriptor in image 1
    :param descriptor 2: Descriptor in image 2
    :param tolerance: Euclidean distance between point in descriptor_1 and point in descriptor_2.

    :returns True: If euclidean distance between each pair of descriptor points is under tolerance.
    """
    for point_order in permutations([0, 1, 2, 3]):
        descriptor_2_perm = np.array(descriptor_2)[list(point_order)]
        distances = np.linalg.norm(descriptor_2_perm - descriptor_1, axis=1)
        if np.all(distances <= tolerance):
            return True
    return False


def match_descriptors(
    descriptors_1, side_lengths_1, descriptors_2, side_lengths_2, tolerance
):
    """
    For each tetrahedral descriptor in reference image, compare it to every descriptor in the target image and pair it with the best match.

    The best match is determined by the difference in the absolute sum of the normalised length of each descriptors's sides, given by get_abs_difference().

    :param descriptors_1: list of descriptors containing their coordinate points, received from get_descriptors().
    :param descriptors_2: list of descriptors containing their coordinate points, received from get_descriptors().
    :param side_lengths_1: normalised side lengths for each descriptor in descriptors_1, received from get_descriptors().
    :param side_lengths_2: normalised side lengths for each descriptor in descriptors_2, received from get_descriptors().
    :param tolerance: Euclidean distance between points in respective image.
    """
    mutual_matches = []
    lowest_differences = []
    difference_matrix = np.zeros(shape=(len(descriptors_1), len(descriptors_2)))

    for i, descriptor_1 in tqdm(
        enumerate(descriptors_1),
        total=len(descriptors_1),
        desc="Matching descriptors - generating the difference matrix",
    ):

        for j, descriptor_2 in enumerate(descriptors_2):
            difference = get_abs_difference(side_lengths_1[i], side_lengths_2[j])
            difference_matrix[i, j] = difference
    # Mutual best match confirmation via difference matrix
    # Need to check that for a given descriptor_1 with best match descriptor_2, that descriptor_2's best match is also descriptor_1
    print("Difference matrix dimensions: ", difference_matrix.shape)
    best_match_1 = np.argmin(difference_matrix, axis=1)
    best_match_2 = np.argmin(difference_matrix, axis=0)

    for i, j in enumerate(best_match_1):
        if best_match_2[j] == i and descriptors_close(
            descriptors_1[i], descriptors_2[j], tolerance
        ):
            mutual_matches.append([descriptors_1[i], descriptors_2[j]])
            lowest_differences.append(difference_matrix[i, j])

    lowest_differences = np.array(lowest_differences)
    sorting_indices = np.argsort(lowest_differences)
    sorted_lowest_differences = lowest_differences[sorting_indices]
    print("Number of mutual matches: ", len(mutual_matches))
    sorted_matches = np.array(mutual_matches, dtype=np.float64)[sorting_indices]

    return sorted_matches, sorted_lowest_differences


def match_descriptor_points(descriptor_1, descriptor_2):
    """
    Matches the vertice coordinates for a pair of matched descriptors.

    Uses the total difference between all possible pairs of matching vertices (6 configurations) as an error metric.

    By selecting the 3 matched point pairs with the lowest error, the correspondence pairs are identified.
    """

    perms = permutations([0, 1, 2, 3])
    lowest_difference = np.inf
    best_permutation = None
    for permutation in perms:
        difference = (
            get_abs_difference(descriptor_1[0], descriptor_2[permutation[0]])
            + get_abs_difference(descriptor_1[1], descriptor_2[permutation[1]])
            + get_abs_difference(descriptor_1[2], descriptor_2[permutation[2]])
            + get_abs_difference(descriptor_1[3], descriptor_2[permutation[3]])
        )
        if difference < lowest_difference:
            lowest_difference = difference
            best_permutation = permutation
    match_1 = [descriptor_1[0], descriptor_2[best_permutation[0]]]
    match_2 = [descriptor_1[1], descriptor_2[best_permutation[1]]]
    match_3 = [descriptor_1[2], descriptor_2[best_permutation[2]]]
    match_4 = [descriptor_1[3], descriptor_2[best_permutation[3]]]
    return match_1, match_2, match_3, match_4


def point_in_pairs_list(point, pairs_list):
    """
    Returns True if a coordinate point exists in a list of coordinate points.

    Used by get_correspondence_pairs() to check if a point has already been paired up (a point will exist in multiple triangles).
    """
    return any(np.array_equal(point, p) for p in pairs_list)


def get_correspondence_pairs(sorted_matched_descriptors):
    """
    Method to get the corresponding point pairs from given matched descriptors.

    Uses the match_descriptor_points() to calculate the correspondence pairs in a matched descriptors pair.

    Repeated for all descriptors, only adding the correspondence pair to the final list if neither point in the correspondence pair has been matched with an earlier point.

    To ensure the best correspondence points are kept, this function should be passed a sorted list of matching descriptors, ordered from best matching to worst matching.
    This way, the earliest correspondence pair contains the points from the best matched descriptors, and we can disregard other possible correspondence pairs.

    :param sorted_matched_descriptors: A list of sorted matched descriptors pairs in order of decreasing similarity, generated from match_descriptors().
    """
    correspondence_pairs = []
    # Lists to keep track of which points have already been matched and are in correspondence_pairs
    used_a = []
    used_b = []
    for descriptor_pair in tqdm(
        sorted_matched_descriptors, desc="Identifying correspondence pairs"
    ):
        match_1, match_2, match_3, match_4 = match_descriptor_points(
            descriptor_pair[0], descriptor_pair[1]
        )
        # Only append correspondence pair to list if both parts of matched pair have not been added as part of another matched pair
        if not point_in_pairs_list(match_1[0], used_a) and not point_in_pairs_list(
            match_1[1], used_b
        ):
            correspondence_pairs.append(match_1)
            used_a.append(match_1[0])
            used_b.append(match_1[1])
        if not point_in_pairs_list(match_2[0], used_a) and not point_in_pairs_list(
            match_2[1], used_b
        ):
            correspondence_pairs.append(match_2)
            used_a.append(match_2[0])
            used_b.append(match_2[1])
        if not point_in_pairs_list(match_3[0], used_a) and not point_in_pairs_list(
            match_3[1], used_b
        ):
            correspondence_pairs.append(match_3)
            used_a.append(match_3[0])
            used_b.append(match_3[1])
        if not point_in_pairs_list(match_4[0], used_a) and not point_in_pairs_list(
            match_4[1], used_b
        ):
            correspondence_pairs.append(match_4)
            used_a.append(match_4[0])
            used_b.append(match_4[1])
    return np.array(correspondence_pairs)
