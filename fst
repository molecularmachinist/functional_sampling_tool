#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from src import cmd_line


if __name__ == '__main__':
    args = cmd_line.argP()
    args.func(args)
