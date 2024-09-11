from nomad.config.models.plugins import ParserEntryPoint
from pydantic import Field

from typing import Optional


class EntryPoint(ParserEntryPoint):
    parser_class_name: str = Field(
        description="""
        The fully qualified name of the Python class that implements the parser.
        This class must have a function `def parse(self, mainfile, archive, logger)`.
    """
    )
    code_name: Optional[str]
    code_homepage: Optional[str]
    code_category: Optional[str]
    metadata: Optional[dict] = Field(
        description="""
        Metadata passed to the UI. Deprecated. """
    )

    def load(self):
        from nomad.parsing import MatchingParserInterface

        return MatchingParserInterface(**self.dict())


h5md_parser_entry_point = EntryPoint(
    name='nomad-parser-h5md',
    aliases=['h5md'],
    description='NOMAD parser for H5MD.',
    python_package='nomad_parser_h5md',
    mainfile_binary_header_re=b'^\\x89HDF',
    mainfile_contents_dict={'__has_all_keys': ['h5md']},
    mainfile_mime_re='(application/x-hdf)',
    mainfile_name_re=r'^.*\.(h5|hdf5)$',
    parser_class_name='nomad_parser_h5md.parsers.parser.H5MDParser',
    code_name='H5MD',
    code_category='MD code',
    metadata={
        'codeCategory': 'MD code',
        'codeLabel': 'H5MD',
        'codeLabelStyle': 'All in capitals',
        'codeName': 'h5md',
        # 'parserDirName': 'dependencies/parsers/atomistic/atomisticparsers/h5md/',
        'parserGitUrl': 'https://github.com/FAIRmat-NFDI/nomad-parser-h5md.git',
        'parserSpecific': '',
        'preamble': '',
        'status': 'production',
        'tableOfFiles': '',
    },
)
