from cif_cli.cli import parser


def test_translate():
    test_parser = parser.parse_args(["transform", "", "", "-t", "0", "0", "0"])
    assert test_parser.translate == [0, 0, 0]


def test_rotate():
    test_parser = parser.parse_args(["transform", "", "", "-r", "0", "0", "0"])
    assert test_parser.rotate == [0, 0, 0]
