# coffeelint: disable=max_line_length

class GroupBrowse
    @controller: (args) ->
        @data = getEndpoint(url="data/meta")
        @tableview = new TableView()
        return @

    @view: (ctrl) ->
        console.log(ctrl.data)
        m(".",
            m.component(TableView, {
                data: ctrl.data
            })
        )