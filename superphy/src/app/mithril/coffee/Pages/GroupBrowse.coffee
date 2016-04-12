# coffeelint: disable=max_line_length

class GroupBrowse extends Page
    @controller: (args) ->
        @data = getEndpoint(url="data/meta")
        @sort_table = (list, attribute = 'data-sort-by') ->
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
        return @
    state = {pageY: 0, pageHeight: window.innerHeight}
    window.addEventListener("scroll", (e) ->
        state.pageY = Math.max(e.pageY || window.pageYOffset, 0)
        state.pageHeight = window.innerHeight
        m.redraw()
    )
    @view: (ctrl) ->
        pageY = state.pageY
        begin = pageY / 60 | 0
        end = begin + (state.pageHeight /60 | 0 + 10)
        offset = pageY % 60
        super m(".Occlusion", {style: {height: ctrl.data().rows.length * 46 + "px", position: "relative", top: -offset + "px"}}, [
            m("table", {style: {top: state.pageY + "px"}}, [
                m("tr", [
                    for header in ctrl.data().headers
                        m('th[data-sort-by=' + header + ']', ctrl.sort_table(list = ctrl.data().rows) ,[header])
                ])
                for row, x in ctrl.data().rows[begin .. end] #when JSON.stringify(row).search(/Unknown/i) > -1
                    m('tr', [
                        for header in ctrl.data().headers
                            m('td', [row[header]])
                    ])
            ])
        ])