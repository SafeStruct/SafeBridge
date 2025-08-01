import logging
import os

class SafeBridgeLogger:
    """
    Singleton logger class for SafeBridge application.
    """
    _instance = None

    def __new__(cls, logfilename: str = "./safebridgeDB/safebridge.log"):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance._logger = get_logger(logfilename, cls._fmt)
            cls._instance._logfilename = logfilename
        return cls._instance

    _fmt = '%(asctime)s | %(levelname)8s | %(filename)s:%(lineno)2d | %(message)s'
    
    def get_logger(self) -> logging.Logger:
        """
        Returns the logger instance.
        """
        return self._logger

    def rename_logfile(self, new_logfilename: str):
        """
        Renames the existing log file to a new file name and updates the logger to use the new file.
        """
        assert isinstance(new_logfilename, str), "new_logfilename must be a string."
        assert new_logfilename, "new_logfilename must not be empty."

        if os.path.exists(self._logfilename):
            os.rename(self._logfilename, new_logfilename)

        self._assing_filename(new_logfilename)

    def _assing_filename(self, logfilename: str):
        """
        Assigns a new log filename to the logger.
        """
        for handler in self._logger.handlers:
            if isinstance(handler, logging.FileHandler):
                handler.close()
                self._logger.removeHandler(handler)

                # Replace with new handler
                new_handler = logging.FileHandler(logfilename, mode='a')
                new_handler.setFormatter(logging.Formatter(self._fmt))
                self._logger.addHandler(new_handler)

                self._logfilename = logfilename
                break

    def existing_logfile(self, existing_logfilename: str):
        """
        Updates the logger to continue writing to an existing log file.
        """
        assert isinstance(existing_logfilename, str), "existing_logfilename must be a string."
        assert existing_logfilename, "existing_logfilename must not be empty."

        self._assing_filename(existing_logfilename)


    
def get_logger(logfilename : str = None, fmt:str = None) -> logging.Logger:
    """
    Initializes and returns a logger instance for the SafeBridge application.
    This function sets up a logger that logs messages to both the console and a specified log file.
    Parameters
    ----------
    logfilename : str
        The name of the log file where logs will be written. Defaults to "safebridge.log".
    Returns
    -------
    logging.Logger
        A configured logger instance that logs messages to both the console and a file.
    Raises
    ------
    AssertionError
        If logfilename is not a string or is None.
    """
    assert(isinstance(logfilename, str)), "logfilename must be a string."
    assert(logfilename is not None), "logfilename must be provided."

    logger = logging.getLogger('safebridgelogger')
    logger.setLevel(logging.DEBUG)

    stdout_handler = get_stream_handler(fmt)
    file_handler = get_file_handler(logfilename, fmt)

    # Add handlers to the logger
    if not logger.hasHandlers():
        logger.addHandler(stdout_handler)
        logger.addHandler(file_handler)
    
    return logger

def get_stream_handler(fmt:str) -> logging.StreamHandler:
    """ Creates and returns a StreamHandler for logging to the console with a custom formatter.
    
    Parameters
    ----------
    fmt : str
        The format string for the log messages.
    Returns
    -------
    logging.StreamHandler
        A configured StreamHandler that logs messages to the console.
    """

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.DEBUG)
    stream_handler.setFormatter(CustomFormatter(fmt))
    return stream_handler

def get_file_handler(logfilename:str, fmt:str) -> logging.FileHandler:
    """ Creates and returns a FileHandler for logging to a file with a custom formatter.

    Parameters
    ----------
    logfilename : str
        The name of the log file where logs will be written.
    fmt : str
        The format string for the log messages.
    Returns
    -------
    logging.FileHandler
        A configured FileHandler that logs messages to the specified file.
    Raises
    ------
    AssertionError
        If logfilename is not a string or is None.
    """
    
    file_handler = logging.FileHandler(filename=logfilename, mode='a')
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(logging.Formatter(fmt))
    return file_handler

class CustomFormatter(logging.Formatter):
    """Logging colored formatter, adapted from https://stackoverflow.com/a/56944256/3638629"""

    grey = '\x1b[38;21m'
    blue = '\x1b[38;5;39m'
    yellow = '\x1b[38;5;226m'
    red = '\x1b[38;5;196m'
    bold_red = '\x1b[31;1m'
    reset = '\x1b[0m'

    def __init__(self, fmt):
        super().__init__()
        self.fmt = fmt
        self.FORMATS = {
            logging.DEBUG: self.grey + self.fmt + self.reset,
            logging.INFO: self.blue + self.fmt + self.reset,
            logging.WARNING: self.yellow + self.fmt + self.reset,
            logging.ERROR: self.red + self.fmt + self.reset,
            logging.CRITICAL: self.bold_red + self.fmt + self.reset
        }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)