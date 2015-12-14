###
    Events
###

events = {}
#Stuctured for Blazegraph response data
events.sort_table = (list, attribute = 'data-sort-by') ->
    { onclick: (e) ->     
        numeric = true
        item = e.target.getAttribute(attribute)
        if item
            first = list[0]
            list.sort (a, b) ->
                if isNaN(parseFloat(a[item]["value"] * 1))
                    if isNaN(parseFloat(b[item]["value"] * 1))
                        if a[item]["value"] > b[item]["value"] then 1 else if b[item]["value"] > a[item]["value"] then -1 else 0
                    else -1
                else if isNaN(parseFloat(b[item]["value"] * 1)) then 1
                else if a[item]["value"] * 1 < b[item]["value"] * 1 then 1 else if b[item]["value"] * 1 < a[item]["value"] * 1 then -1 else 0

            if first == list[0]
                list.reverse()
        return
    }

m.route(document.body, "/", {
    "/": Home.get()
    "/home": Foo.get()
    "/meta": GroupBrowse.get()
    "/gbrowse": GroupBrowse.get()
    "/groups": GroupBrowse.get()
})