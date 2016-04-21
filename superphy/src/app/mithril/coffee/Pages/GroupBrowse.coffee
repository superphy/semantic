# coffeelint: disable=max_line_length

class GroupBrowse extends Page
    Routes.add('/gbrowse', @)

    @controller: (args) ->
        @data = getEndpoint(url="data/meta")()
        return @
        
    @view: (ctrl) ->
        return super(
            m(".",
                m.component(TableView, {
                    data: ctrl.data
                })
            )
        )