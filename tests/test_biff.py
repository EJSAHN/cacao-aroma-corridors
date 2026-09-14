from __future__ import annotations
import struct
import unittest
from cacao_marker_harmonization.io.biff8 import BiffBook, decode_sst
from cacao_marker_harmonization.sources.common import DataError

class BiffPrimitiveTests(unittest.TestCase):

    def test_rk_integer(self):
        self.assertEqual(BiffBook.rk(42 << 2 | 2), 42)

    def test_rk_divided(self):
        self.assertEqual(BiffBook.rk(1234 << 2 | 3), 12.34)

    def test_single_byte_sst(self):
        b = struct.pack('<II', 1, 1) + struct.pack('<HB', 5, 0) + b'hello'
        self.assertEqual(decode_sst([b]), ['hello'])

    def test_sst_continuation(self):
        b = struct.pack('<II', 1, 1) + struct.pack('<HB', 5, 0) + b'he'
        self.assertEqual(decode_sst([b, b'\x00llo']), ['hello'])

    def test_sst_wide_string(self):
        word = 'fruité'
        b = struct.pack('<II', 1, 1) + struct.pack('<HB', len(word), 1) + word.encode('utf-16le')
        self.assertEqual(decode_sst([b]), [word])

    def test_truncated_sst(self):
        b = struct.pack('<II', 1, 1) + struct.pack('<HB', 5, 0) + b'he'
        with self.assertRaises(DataError):
            decode_sst([b])
