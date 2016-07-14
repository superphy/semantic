# FOR Mockcomponent

# coffeelint: disable=max_line_length

class Mock extends Page
    Routes.add("?/mock", @)
   
    @controller: (args) ->
        @data = getEndpoint2(url="data/meta")
        return @
    @view: (ctrl) ->
        super m("."
            m.component(Mockcomponent, data: ctrl.data)
        )
