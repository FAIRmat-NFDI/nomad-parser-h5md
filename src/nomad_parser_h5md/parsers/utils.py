import functools
import os
import re
from collections.abc import Callable
from glob import glob
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from structlog.stdlib import (
        BoundLogger,
    )
from nomad.metainfo import Section, SubSection
from nomad.utils import get_logger


def remove_mapping_annotations(property: Section, max_depth: int = 5) -> None:
    """
    Remove mapping annotations from the input section definition, all its quantities
    and sub-sections recursively.

    Args:
        property (Section): The section definition to remove the annotations from.
        max_depth (int, optional): The maximum depth of the recursion for sub-sections
            using the same section as parent.
    """

    def _remove(property: Section | SubSection, depth: int = 0):
        if depth > max_depth:
            return

        annotation_key = 'mapping'
        property.m_annotations.pop(annotation_key, None)

        depth += 1
        property_section = (
            property.sub_section if isinstance(property, SubSection) else property
        )
        for quantity in property_section.all_quantities.values():
            quantity.m_annotations.pop(annotation_key, None)

        for sub_section in property_section.all_sub_sections.values():
            if sub_section.m_annotations.get(annotation_key):
                _remove(sub_section, depth)
            elif sub_section.sub_section.m_annotations.get(annotation_key):
                _remove(sub_section.sub_section, depth)
            else:
                for (
                    inheriting_section
                ) in sub_section.sub_section.all_inheriting_sections:
                    if inheriting_section.m_annotations.get(annotation_key):
                        _remove(inheriting_section, depth)

    _remove(property)
