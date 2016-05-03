# coffeelint: disable=max_line_length

class GroupBrowse extends Page
    Routes.add('/gbrowse', @)

    @controller: (args) ->
        @data = getEndpoint2(url="data/meta")()
        return @
        
    @view: (ctrl) ->
        super m("."
            m.component(Table, data: ctrl.data)
        )
