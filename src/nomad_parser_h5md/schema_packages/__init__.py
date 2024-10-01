from nomad.config.models.plugins import SchemaPackageEntryPoint
from pydantic import Field


class H5MDSchemaPackageEntryPoint(SchemaPackageEntryPoint):
    parameter: int = Field(0, description='Custom configuration parameter')

    def load(self):
        from nomad_parser_h5md.schema_packages.schema_package import m_package

        return m_package


h5md_schema_package_entry_point = H5MDSchemaPackageEntryPoint(
    name='H5MDSchemaPackage',
    description='H5MD schema package entry point configuration.',
)
