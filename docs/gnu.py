#!/usr/bin/env python
# -*- coding: utf-8 -*-

exec(open("../pydna/_version.py").read())
release = get_versions()["version"]
version = '.'.join(release.split('.')[:2])

print(version)
