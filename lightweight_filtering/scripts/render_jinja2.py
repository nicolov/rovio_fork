#!/usr/bin/env python

from __future__ import print_function

import re
import os
import glob
import argparse

from jinja2 import Environment, FileSystemLoader


def generate(input_dir, output_dir):
    ''' Renders templates that match input_dir/*.ext.j2, producing the
    corresponding output in output_dir/*.ext.

    Uses input_dir as the path for the template loader, and skips input files
    that start with an underscore. '''

    template_loader_dirs = [
        os.path.join(os.path.dirname(__file__), '../templates'),
        input_dir,
    ]

    env = Environment(
        loader=FileSystemLoader(template_loader_dirs),
        trim_blocks=True,
        lstrip_blocks=True)

    # Build list of template paths
    template_paths = glob.glob(os.path.join(input_dir, '*.j2'))
    should_skip = lambda path: not os.path.basename(path).startswith('_')
    template_paths = filter(should_skip, template_paths)

    for tpl_path in template_paths:
        tpl_name = os.path.basename(tpl_path)
        # Strip .j2 if needed
        output_name = re.sub(r'\.j2$', '', tpl_name)
        output_path = os.path.join(output_dir, output_name)

        template = env.get_template(tpl_name)

        with open(output_path, 'w') as fp:
            fp.write(template.render())


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-dir', help='Path to input directory', required=True)
    parser.add_argument('-o', '--output-dir', help='Path to output directory', required=True)
    args = parser.parse_args()

    generate(args.input_dir, args.output_dir)
