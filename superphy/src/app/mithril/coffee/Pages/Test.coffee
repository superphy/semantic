# FOR TABLE 2 BELOW. TO GET DATA FROM DATABASE, USE TEMPLATE BELOW

# coffeelint: disable=max_line_length

class Test extends Page
    Routes.add("/test", @)
   
    @controller: (args) ->
        @data = getEndpoint2(url="data/meta")
        return @

    #To be removed as soon as offset on page is working
    @view: (ctrl) ->
        super m("."
            m.component(Table_2, data: ctrl.data)
        )
