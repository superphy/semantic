# coffeelint: disable=max_line_length

class GroupBrowse
    @controller: (args) ->
        @data = getEndpoint(url="data/meta")()
        return @

    @view: (ctrl) ->
        m(".",
            m.component(TableView, {
                data: ctrl.data
            })
        )