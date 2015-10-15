from sqlalchemy import *
from migrate import *


from migrate.changeset import schema
pre_meta = MetaData()
post_meta = MetaData()
genome = Table('genome', post_meta,
    Column('id', Integer, primary_key=True, nullable=False),
)

user = Table('user', post_meta,
    Column('id', Integer, primary_key=True, nullable=False),
    Column('email', String(length=254)),
    Column('password', String(length=128)),
    Column('first_name', String(length=255)),
    Column('last_name', String(length=255)),
    Column('org', String(length=255)),
)


def upgrade(migrate_engine):
    # Upgrade operations go here. Don't create your own engine; bind
    # migrate_engine to your metadata
    pre_meta.bind = migrate_engine
    post_meta.bind = migrate_engine
    post_meta.tables['genome'].create()
    post_meta.tables['user'].create()


def downgrade(migrate_engine):
    # Operations to reverse the above upgrade go here.
    pre_meta.bind = migrate_engine
    post_meta.bind = migrate_engine
    post_meta.tables['genome'].drop()
    post_meta.tables['user'].drop()
