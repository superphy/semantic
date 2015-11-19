#!/usr/bin/env python

"""mapper

Converts input metadata JSON input with varisous metadata type-value pairs
to URIs used in the Superphy DB

Example:
    Examples can be given using either the ``Example`` or ``Examples``
    sections. Sections support any reStructuredText formatting, including
    literal blocks::

        $ python example_google.py

Section breaks are created by resuming unindented text. Section breaks
are also implicitly created anytime a new section starts.


"""

import logging
from superphy.config import parser

__author__ = "Matt Whiteside"
__copyright__ = "Copyright 2015, Public Health Agency of Canada"
__license__ = "CC-SY"
__version__ = "4.0"
__maintainer__ = "Matt Whiteside"
__email__ = "matthew.whiteside@canada.ca"


class SuperphyMapperError(Exception):
    """Exceptions are documented in the same way as classes.

    """

    def __init__(Exception):
        pass


class Mapper(object):
    """Maps various metadata JSON inputs to SuperphyDB URIs

    If the class has public attributes, they may be documented here
    in an ``Attributes`` section and follow the same formatting as a
    function's ``Args`` section. Alternatively, attributes may be documented
    inline with the attribute's declaration (see __init__ method below).

    Properties created with the ``@property`` decorator should be documented
    in the property's getter method.

    Attribute and property types -- if given -- should be specified according
    to `PEP 484`_, though `PEP 484`_ conformance isn't required or enforced.

    Attributes:
        attr1 (str): Description of `attr1`.
        attr2 (Optional[int]): Description of `attr2`.


    .. _PEP 484:
       https://www.python.org/dev/peps/pep-0484/

    """

    def __init__(self, logger=None):
        """


        Args:
            param1 (str): Description of `param1`.
            param2 (Optional[int]): Description of `param2`. Multiple
                lines are supported.
            param3 (List[str]): Description of `param3`.

        """

        self.logger = logger or logging.getLogger(__name__)

        # Retrieve superphy config
        self.config = parser.read()

        # Init Rdflib endpoint 