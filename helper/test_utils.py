import tempfile
import unittest
from pathlib import Path

from helper.utils import count_lines


class UtilsTests(unittest.TestCase):
    """Test utils"""

    def test_count_lines(self):
        """Test count lines"""

        with tempfile.NamedTemporaryFile(mode="wt") as file:
            path = Path(file.name)
            file.write("line1\nline2\nline3\n\nline5\n")
            file.flush()
            self.assertEqual(count_lines(path), 5)
            file.truncate(0)
            file.seek(0)
            self.assertEqual(count_lines(path), 0)
