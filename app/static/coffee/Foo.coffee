class List
    model: () =>
        @strings= ['Ecoli']
        @value= m.prop('')
        return
    constructor: ->
        @model()
        return
    view: () =>
        return [
            #List of objects in the array
            m('ol', {}, [
                m('li',[
                    m('button[index=' + index + ']',events.delete_item(@strings),"DELETE")
                    m('text',string)
                ]) for string, index in @strings
            ])
            #Form to input objects in the array
            m('label', ['New Thing: '
                m('input', { oninput: m.withAttr('value', @value), value : @value() })
                m('button', { disabled : !@value(), onclick: => if @value() then @strings.push(@value()) and @value("")},["Add"])
            ])
        ]

class Searchbar
    model: () =>
    controller: () ->

    view: () =>
        return [
            #m('th[data-sort-by=' + header + ']',events.sort_table(list = data.genomes) ,[header]) 
            m('input[data-search]',{},[])
        ]
class Bar
    data = {}
    data.headers= ["A", "B", "C", "D"]
    data.genomes= [
                    {"A": {"value": "AB"}, "B": {"value": "15572"}, "C": {"value":""},"D": {"value": "FOO"}},
                    {"A": {"value": "ABC"}, "B": {"value": "15572"}, "C": {"value":"2435852"},"D": {"value": "TEST"}},
                    {"A": {"value": "ABCD"}, "B": {"value": "15572"}, "C": {"value":"HELLO"},"D": {"value": "WORLD"}},
                    {"A": {"value": "ABCDEF"}, "B": {"value": "15572"}, "C": {"value":"2435852"},"D": {"value": "BAR"}},
                    {"A": {"value": "ABCDEF"}, "B": {"value": "15572"}, "C": {"value":"2435852"},"D": {"value": "BAR"}},
                    {"A": {"value": "ABCDEF"}, "B": {"value": "15572"}, "C": {"value":"2435852"},"D": {"value": "BAR"}},
                    {"A": {"value": "ABCDEF"}, "B": {"value": "15572"}, "C": {"value":"2435852"},"D": {"value": "BAR"}},
                    {"A": {"value": "ABCDEF"}, "B": {"value": "15572"}, "C": {"value":"2435852"},"D": {"value": "BAR"}},
                    {"A": {"value": "ABCDEF"}, "B": {"value": "15572"}, "C": {"value":"2435852"},"D": {"value": "BAR"}},
                    {"A": {"value": "ABCDEF"}, "B": {"value": "15572"}, "C": {"value":"2435852"},"D": {"value": "BAR"}},
                    {"A": {"value": "ABCDEF"}, "B": {"value": "15572"}, "C": {"value":"2435852"},"D": {"value": "BAR"}},
                    {"A": {"value": "ABCDEF"}, "B": {"value": "15572"}, "C": {"value":"2435852"},"D": {"value": "BAR"}},
                    {"A": {"value": "ABCDEF"}, "B": {"value": "15572"}, "C": {"value":"2435852"},"D": {"value": "BAR"}}]
    response: () =>
        return data


class Foo extends Page_Template
    model: () =>
        @search ?= new Searchbar()
        @data ?= new Bar
        @table ?= new Table(@data)
        return
    controller: (options) =>
        @model()
        return
    view: () =>
        return [
            @search.view()
            Header.getView()
            @table.view()
        ]