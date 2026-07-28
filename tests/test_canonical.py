from nfc_cryst.canonical import canonical_bytes, semantic_sha256


def test_canonical_json_is_order_independent() -> None:
    first = {"b": 2, "a": [1, 3]}
    second = {"a": [1, 3], "b": 2}
    assert canonical_bytes(first) == canonical_bytes(second)
    assert semantic_sha256(first) == semantic_sha256(second)
