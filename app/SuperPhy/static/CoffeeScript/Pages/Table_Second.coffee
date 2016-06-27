#Started here
class Foo_ooo extends Page
    Routes.add('/foo', @)
    @controller: (args) ->
        return @
    @view: (ctrl, args) ->
        super m.component(Table_second)

# End
class Table_second