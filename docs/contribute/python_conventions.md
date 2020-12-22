# Python conventions and style guide



## typing

all python code must have type hints. 

#### PEP-526:
```
from typing import List, Set, Dict, Tuple, Optional

x: int = 1
x: float = 1.0
x: bool = True
x: str = "test"
x: bytes = b"test"

x: List[int] = [1]
x: Set[int] = {6, 7}

x: Dict[str, float] = {'field': 2.0}

x: Tuple[int, str, float] = (3, "yes", 7.5)

# For tuples of variable size, we use one type and ellipsis
x: Tuple[int, ...] = (1, 2, 3)

# Use Optional[] for values that could be None
x: Optional[str] = some_function()

# functions
def f(num1: int, num2: int, my_float: float = 3.5) -> float:
    return num1 + num2 + my_float

# callable (function) value
x: Callable[[int, float], float] = f

# function that returns iterators
def returns_iterator(n: int) -> Iterator[int]:
    i = 0
    while i < n:
        yield i
        i += 1

 def returns_none(__x: int) -> None:
    pass


from typing import Union, Any, List, Optional, cast

# Use Union when something could be one of a few types
x: List[Union[int, str]] = [3, 5, "test", "fun"]


# Use Any if you don't know the type of something or it's too
# dynamic to write a type for
x: Any = mystery_function()


# This makes each positional arg and each keyword arg a "str"
def call(self, *args: str, **kwargs: str) -> str:
    request = make_request(*args, **kwargs)
    return self.do_api_query(request)

# type: ignore comment to skip
x = confusing_function()  # type: ignore

```

## naming:
- function: lowercase word(s). underscore separated
- variable: lowercase word(s). underscore separated
- class: CamelCase (camelCase is not acceptable)
- constant: UPPERCASE
- module: lowercase word(s). underscore separated

## blank lines:
- surround top level functions with 2 lines
- surround method definitions inside classes with a single blank line
- use blank line in function to separate code when necessary

## Indent:
4 spaces. no exceptions.


## Comments:

- always add a space after `#`
```
	for i in range(0, 10):
    # Loop over i ten times and print out the value of i, followed by a
    # new line character
        print(i, '\n')
    x: int = 5  # This is an inline comment
```
## Strings

- all strings should use `"`quotes
- always use fstrings
- when using fstrings, key access should use `'`

Example:
```
f"here is {mydict['key']}"
```


## line breaks

- line length under 80
- use hanging indents.
```
age: Dict[str, int] = {
    'Alice': 24,
    'Bob': 28,
    'Ann': 26,
}
def function(
        arg_one: int,
        arg_two: int,
        arg_three: int,
        arg_four: int):
    return arg_one
```

- operators on the new line
```
total = (first_variable
         + second_variable
         - third_variable)
```



## imports

yes:
```
from mypkg import example1, \
    example2, example3
```

```
from mypkg import example1
from mypkg import example2
from mypkg import example3
from mypkg import example4

```

no:
```
import sys, os
```
```
from sys import *
```

## Doc strings

are required. 

```
def get_read_tag(fastq_read: str):
    """
    explain function here  
    Parameters
    ----------
    fastq_read: str
        file path


    Returns
    -------
    tag: str
        read tag
    """
```


**Class docstrings:**

The docstring for a class should summarize its behavior and list the public methods and instance variables. If the class is intended to be subclassed, and has an additional interface for subclasses, this interface should be listed separately (in the docstring). The class constructor should be documented in the docstring for its __init__ method. Individual methods should be documented by their own docstring.

```
class Person:
    """
    A class to represent a person.

    ...

    Attributes
    ----------
    name : str
        first name of the person
    surname : str
        family name of the person
    age : int
        age of the person

    Methods
    -------
    info(additional=""):
        Prints the person's name and age.
    """

    def __init__(self, name: str, surname: str, age: int):
        """
        Constructs all the necessary attributes for the person object.

        Parameters
        ----------
            name : str
                first name of the person
            surname : str
                family name of the person
            age : int
                age of the person
        """

        self.name = name
        self.surname = surname
        self.age = age

    def info(self, additional: str=""):
        """
        Prints the person's name and age.
        If the argument 'additional' is passed, then it is appended after the main info.

        Parameters
        ----------
        additional : str, optional
            More info to be displayed (default is None)

        Returns
        -------
        None
        """

        print(f'My name is {self.name} {self.surname}. I am {self.age} years old.' + additional)
```

## Guides
- always follow [PEP-8](https://www.python.org/dev/peps/pep-0008/)
- refer to [mypy](https://mypy.readthedocs.io/en/stable/cheat_sheet_py3.html) for type hints

