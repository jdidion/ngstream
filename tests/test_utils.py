from ngstream.utils import *


def test_index_batcher():
    batcher = IndexBatcher(batch_size=10)
    assert list(batcher(100)) == [
        (i, start, 10) for i, start in enumerate(range(0, 100, 10))
    ]

    batcher = IndexBatcher(batch_start=1, batch_size=10, batch_step=4)
    assert list(batcher(100)) == [(0, 10, 10), (1, 50, 10), (2, 90, 10)]

    batcher = IndexBatcher(
        item_start=5, item_stop=95, batch_start=1, batch_size=10, batch_step=4
    )
    assert list(batcher(100)) == [(0, 15, 10), (1, 55, 10)]

    batcher = IndexBatcher(
        item_start=5,
        item_stop=95,
        item_limit=15,
        batch_start=1,
        batch_size=10,
        batch_step=4,
    )
    assert list(batcher(100)) == [(0, 15, 10), (1, 55, 5)]
