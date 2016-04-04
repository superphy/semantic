###
    Events
###

events = {}
events.sort_table = (list, attribute = 'data-sort-by') ->
    { onclick: (e) ->
        item = e.target.getAttribute(attribute)
        if item
            first = list[0]
            list.sort (a, b) ->
                if isNaN(parseFloat(a[item] * 1))
                    if isNaN(parseFloat(b[item] * 1))
                        if a[item] > b[item] then 1
                        else if b[item] > a[item] then -1
                        else 0
                    else -1
                else if isNaN(parseFloat(b[item] * 1)) then 1
                else if a[item] * 1 < b[item] * 1 then 1
                else if b[item] * 1 < a[item] * 1 then -1
                else 0
            if first == list[0]
                list.reverse()
        return
    }
events.delete_item = (list) ->
    {onclick: (e) ->
        console.log('*****************')
        index = e.target.getAttribute('index')
        for x,y in list
            console.log("x: " + x, "y: " + y, "index: " + index)
            if parseInt(index) == parseInt(y)
                list.splice(y,1)

        console.log("list: " + list)
    }

m.route(document.body, "/", {
    "/": {view: -> [
        m.component(Navbar)
        m.component(Home)
    ]}
    "/home": {view: -> [
        m.component(Navbar)
    ]}
    "/meta": {view: -> [
        m.component(Navbar)
        m.component(GroupBrowseFirst10)
    ]}
    "/gbrowse": {view: -> [
        m.component(Navbar)
        m.component(GroupBrowse)
    ]}
    "/groups": {view: -> [
        m.component(Navbar)
        m.component(ComponentTemplate)
    ]}
    "/factors": {view: -> [
        m.component(Navbar)
        m.component(FooFoo, {id: 1})
        m.component(FooFoo, {id: 2})
        m.component(FooFoo, {id: 3})
    ]}
})