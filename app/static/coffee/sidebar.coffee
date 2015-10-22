#model
Thing =
    list: -> ['one','two','three']

#top level component
sidebar =
    controller: ->
        @list = new list.controller({
            visible: ->
                item.name.indexOf(@filter.searchTerm()) > -1
        })
        @filter = new filter.controller()
    view: (ctrl) ->
        m(".row", [
            m(".col-md-2", [
                filter.view(ctrl.filter)
            ])
            m(".col-md-10", [
                list.view(ctrl.list)
            ])
        ])

#filter
filter =
    controller: ->
        @searchTerm = m.prop("")
    view: ->
        m("input", {oninput: m.withAttr("value", ctrl.searchTerm)})
        
list =
    controller: ->
        @items = Thing.list()
        @visible = options.visible

    view: (ctrl) ->
        m("table", [
            m("tr", [
                m("td", item.id)
                m("td", item.name)]) for item in(
                    ctrl.items().filter(ctrl.visible))
            ])





