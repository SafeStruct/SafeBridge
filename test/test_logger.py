import unittest
import os
import tempfile
import logging
from safebridge.logger import SafeBridgeLogger, get_logger


class TestSafeBridgeLogger(unittest.TestCase):
    
    def setUp(self):
        """Set up test fixtures before each test method."""
        self.temp_dir = tempfile.mkdtemp()
        self.test_logfile = os.path.join(self.temp_dir, "test.log")
        self.test_format = "%(levelname)s - %(message)s"
        
    def tearDown(self):
        """Clean up after each test method."""
        # Clean up any created log files
        for file in os.listdir(self.temp_dir):
            try:
                os.remove(os.path.join(self.temp_dir, file))
            except:
                pass
        try:
            os.rmdir(self.temp_dir)
        except:
            pass
        
        # Clear all handlers from logger
        logger = logging.getLogger('SafeBridge')
        for handler in logger.handlers[:]:
            handler.close()
            logger.removeHandler(handler)

    def test_init_default_parameters(self):
        """Test SafeBridgeLogger initialization with default parameters."""
        logger_instance = SafeBridgeLogger()
        
        self.assertEqual(logger_instance._logfilename, "safebridge.log")
        self.assertEqual(logger_instance._fmt, '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        self.assertIsInstance(logger_instance._logger, logging.Logger)
        self.assertEqual(logger_instance._logger.name, 'SafeBridge')
        self.assertEqual(logger_instance._logger.level, logging.DEBUG)

    def test_init_custom_parameters(self):
        """Test SafeBridgeLogger initialization with custom parameters."""
        logger_instance = SafeBridgeLogger(self.test_logfile, self.test_format)
        
        self.assertEqual(logger_instance._logfilename, self.test_logfile)
        self.assertEqual(logger_instance._fmt, self.test_format)
        self.assertTrue(os.path.exists(self.test_logfile))

    def test_handlers_setup(self):
        """Test that both file and console handlers are properly set up."""
        logger_instance = SafeBridgeLogger(self.test_logfile)
        handlers = logger_instance._logger.handlers
        
        self.assertEqual(len(handlers), 2)
        
        # Check for FileHandler
        file_handlers = [h for h in handlers if isinstance(h, logging.FileHandler)]
        self.assertEqual(len(file_handlers), 1)
        
        # Check for StreamHandler (console)
        stream_handlers = [h for h in handlers if isinstance(h, logging.StreamHandler) 
                          and not isinstance(h, logging.FileHandler)]
        self.assertEqual(len(stream_handlers), 1)

    def test_get_logger(self):
        """Test get_logger method returns the correct logger instance."""
        logger_instance = SafeBridgeLogger(self.test_logfile)
        returned_logger = logger_instance.get_logger()
        
        self.assertIs(returned_logger, logger_instance._logger)
        self.assertEqual(returned_logger.name, 'SafeBridge')

    def test_rename_logfile_existing_file(self):
        """Test renaming an existing log file."""
        # Create initial logger and write to file
        logger_instance = SafeBridgeLogger(self.test_logfile)
        logger_instance._logger.info("Test message")
        
        # Verify original file exists
        self.assertTrue(os.path.exists(self.test_logfile))
        
        # Rename the file
        new_logfile = os.path.join(self.temp_dir, "renamed.log")
        logger_instance.rename_logfile(new_logfile)
        
        # Verify old file is renamed and new file exists
        self.assertFalse(os.path.exists(self.test_logfile))
        self.assertTrue(os.path.exists(new_logfile))
        self.assertEqual(logger_instance._logfilename, new_logfile)

    def test_rename_logfile_nonexistent_file(self):
        """Test renaming when original file doesn't exist."""
        nonexistent_file = os.path.join(self.temp_dir, "nonexistent.log")
        logger_instance = SafeBridgeLogger(nonexistent_file)
        
        # Remove the file if it was created
        if os.path.exists(nonexistent_file):
            os.remove(nonexistent_file)
        
        new_logfile = os.path.join(self.temp_dir, "new.log")
        logger_instance.rename_logfile(new_logfile)
        
        self.assertEqual(logger_instance._logfilename, new_logfile)

    def test_rename_logfile_invalid_input(self):
        """Test rename_logfile with invalid inputs."""
        logger_instance = SafeBridgeLogger(self.test_logfile)
        
        # Test with non-string input
        with self.assertRaises(AssertionError):
            logger_instance.rename_logfile(123)
        
        # Test with empty string
        with self.assertRaises(AssertionError):
            logger_instance.rename_logfile("")

    def test_existing_logfile_valid(self):
        """Test existing_logfile with a valid existing file."""
        # Create a log file first
        existing_file = os.path.join(self.temp_dir, "existing.log")
        with open(existing_file, 'w') as f:
            f.write("Existing log content\n")
        
        logger_instance = SafeBridgeLogger(self.test_logfile)
        logger_instance.existing_logfile(existing_file)
        
        self.assertEqual(logger_instance._logfilename, existing_file)

    def test_existing_logfile_invalid_input(self):
        """Test existing_logfile with invalid inputs."""
        logger_instance = SafeBridgeLogger(self.test_logfile)
        
        # Test with non-string input
        with self.assertRaises(AssertionError):
            logger_instance.existing_logfile(123)
        
        # Test with empty string
        with self.assertRaises(AssertionError):
            logger_instance.existing_logfile("")
        
        # Test with non-existent file
        nonexistent_file = os.path.join(self.temp_dir, "nonexistent.log")
        with self.assertRaises(AssertionError):
            logger_instance.existing_logfile(nonexistent_file)

    def test_assign_filename(self):
        """Test _assign_filename private method."""
        logger_instance = SafeBridgeLogger(self.test_logfile)
        new_file = os.path.join(self.temp_dir, "assigned.log")
        
        logger_instance._assign_filename(new_file)
        
        self.assertEqual(logger_instance._logfilename, new_file)
        # Verify file handler is updated
        file_handlers = [h for h in logger_instance._logger.handlers 
                        if isinstance(h, logging.FileHandler)]
        self.assertEqual(len(file_handlers), 1)

class TestGetLoggerFunction(unittest.TestCase):
    
    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        
    def tearDown(self):
        """Clean up after each test."""
        # Clean up log files
        for file in os.listdir(self.temp_dir):
            try:
                os.remove(os.path.join(self.temp_dir, file))
            except:
                pass
        try:
            os.rmdir(self.temp_dir)
        except:
            pass
        
        # Clear logger handlers
        logger = logging.getLogger('SafeBridge')
        for handler in logger.handlers[:]:
            handler.close()
            logger.removeHandler(handler)

    def test_get_logger_default_parameters(self):
        """Test get_logger function with default parameters."""
        logger = get_logger()
        
        self.assertIsInstance(logger, logging.Logger)
        self.assertEqual(logger.name, 'SafeBridge')
        self.assertTrue(os.path.exists("safebridge.log"))

    def test_get_logger_custom_parameters(self):
        """Test get_logger function with custom parameters."""
        custom_logfile = os.path.join(self.temp_dir, "custom.log")
        custom_format = "%(levelname)s: %(message)s"
        
        logger = get_logger(custom_logfile, custom_format)
        
        self.assertIsInstance(logger, logging.Logger)
        self.assertEqual(logger.name, 'SafeBridge')
        self.assertTrue(os.path.exists(custom_logfile))

    def test_get_logger_none_logfilename(self):
        """Test get_logger function with None logfilename."""
        logger = get_logger(logfilename=None)
        
        self.assertIsInstance(logger, logging.Logger)
        self.assertTrue(os.path.exists("safebridge.log"))



if __name__ == '__main__':
    unittest.main()