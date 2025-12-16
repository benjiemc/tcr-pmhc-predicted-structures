"""Functions for padding things."""

from functools import lru_cache


@lru_cache
def centre_pad(length: int, pad_length: int) -> tuple[str, int, int]:
    """Determine the amount of padding or cropping to centre something with the specified length."""
    if length < pad_length:
        if (pad_length % 2) == (length % 2):
            pad_per_side = (pad_length - length) // 2
            pad_right_side = pad_per_side
            pad_left_side = pad_per_side

        else:
            pad_right_side = (pad_length - length) // 2
            pad_left_side = pad_right_side + 1

        return 'pad', pad_left_side, pad_right_side

    if length > pad_length:
        if (pad_length % 2) == (length % 2):
            crop_per_side = (length - pad_length) // 2
            crop_right_side = crop_per_side
            crop_left_side = crop_per_side

        else:
            crop_right_side = (length - pad_length) // 2
            crop_left_side = crop_right_side + 1

        return 'crop', crop_left_side, crop_right_side

    return '', 0, 0


@lru_cache
def right_pad(length: int, pad_length: int) -> tuple[str, int, int]:
    """Right pad or crop things based on the specified pad length.

    Examples:
        >>> right_pad(3, 5)
        'pad', 0, 2

        >>> right_pad(3, 2)
        'crop', 0, 1

        >>> right_pad(3, 3)
        '', 0, 0

    """
    if length < pad_length:
        num_to_add = pad_length - length
        return 'pad', 0, num_to_add

    if length > pad_length:
        num_to_remove = length - pad_length
        return 'crop', 0, num_to_remove

    return '', 0, 0


@lru_cache
def left_pad(length: int, pad_length: int) -> tuple[str, int, int]:
    """Left pad or crop things based on the specified pad length.

    Examples:
        >>> left_pad(3, 5)
        'pad', 2, 0

        >>> left_pad(3, 2)
        'crop', 1, 0

        >>> left_pad(['A', 'B', 'C'], 3)
        'none', 0, 0

    """
    if length < pad_length:
        num_to_add = pad_length - length
        return 'pad', num_to_add, 0

    if length > pad_length:
        num_to_remove = length - pad_length
        return 'crop', num_to_remove, 0

    return 'none', 0, 0
