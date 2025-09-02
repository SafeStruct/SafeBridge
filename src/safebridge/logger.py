import os
import logging

class SafeBridgeLogger:
    def __init__(self, logfilename: str = "safebridge.log", fmt: str = None):
        """
        Initialize the SafeBridge logger.
        
        Parameters
        ----------
        logfilename : str
            The name of the log file
        fmt : str
            The logging format string

        Methods
        -------
        rename_logfile(new_logfilename: str):
            Renames the existing log file to a new file name and updates the logger to use the new file.
        existing_logfile(existing_logfilename: str):
            Updates the logger to continue writing to an existing log file.
        get_logger() -> logging.Logger:
            Returns the logger instance.
        
        """
        self._logfilename = logfilename
        self._fmt = fmt or '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        self._logger = logging.getLogger('SafeBridge')
        self._logger.setLevel(logging.DEBUG)
        
        # Clear existing handlers to prevent duplicates
        if self._logger.handlers:
            self._logger.handlers.clear()

        # Add file handler
        file_handler = logging.FileHandler(logfilename, mode='a')
        file_handler.setFormatter(logging.Formatter(self._fmt))
        self._logger.addHandler(file_handler)
        
        # Add console handler
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(logging.Formatter(self._fmt))
        self._logger.addHandler(console_handler)
    
    def rename_logfile(self, new_logfilename: str):
        """ Renames the existing log file to a new file name and updates the logger to use the new file.

        Parameters
        ----------
        new_logfilename : str
            The new name for the log file.

        Raises
        ------
        AssertionError
            If new_logfilename is not a string or is empty.
        """
        assert isinstance(new_logfilename, str), "new_logfilename must be a string."
        assert new_logfilename, "new_logfilename must not be empty."

        if os.path.exists(self._logfilename):
            os.rename(self._logfilename, new_logfilename)

        self._assign_filename(new_logfilename)

    def _assign_filename(self, logfilename: str):
        """ Assigns a new log filename to the logger.

        Parameters
        ----------
        logfilename : str
            The new log filename to assign.

        Raises
        ------
        AssertionError
            If logfilename is not a string or is empty and does not exist.
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
        """ Updates the logger to continue writing to an existing log file.

        Parameters
        ----------
        existing_logfilename : str
            The name of the existing log file to continue writing to.
        
        Raises
        ------
        AssertionError
            If existing_logfilename is not a string, is empty, or does not exist.
        """
        assert isinstance(existing_logfilename, str), "existing_logfilename must be a string."
        assert existing_logfilename, "existing_logfilename must not be empty."
        assert os.path.exists(existing_logfilename), f"Log file {existing_logfilename} does not exist."

        self._assign_filename(existing_logfilename)

    def get_logger(self) -> logging.Logger:
        """
        Returns the logger instance.
        """
        return self._logger

def get_logger(logfilename: str = None, fmt: str = None) -> logging.Logger:
    """
    Initializes and returns a logger instance for the SafeBridge application.
    This function sets up a logger that logs messages to both the console and a specified log file.
    
    Parameters
    ----------
    logfilename : str, optional
        The name of the log file where logs will be written. Defaults to "safebridge.log".
    fmt : str, optional
        The logging format string. Uses default format if None.
        
    Returns
    -------
    logging.Logger
        A configured logger instance that logs messages to both the console and a file.
    """
    if logfilename is None:
        logfilename = "safebridge.log"
    
    bridge_logger = SafeBridgeLogger(logfilename, fmt)
    return bridge_logger.get_logger()