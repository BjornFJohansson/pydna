#! /usr/bin/python
import logging
import test_log


def main():
    logging.error("hello")


if __name__ == "__main__":
    logging.basicConfig(filename="example.log", filemode="a", level=logging.DEBUG)

    main()
