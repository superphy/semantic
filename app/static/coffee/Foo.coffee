class Searchbar
    model: () =>
    controller: () ->

    view: () =>
        return [
            #m('th[data-sort-by=' + header + ']',events.sort_table(list = data.genomes) ,[header]) 
            m('input[data-search',{},[])
        ]
class Bar
    response: () ->
        return {
            "headers": ["A", "B", "C", "D"],"genomes": [
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
                {"A": {"value": "ABCDEF"}, "B": {"value": "15572"}, "C": {"value":"2435852"},"D": {"value": "BAR"}}]}

class Foo extends Page_Template
    model: () =>
        @data ?= new Bar
        @table ?= new Table(@data)
        return
    controller: (options) =>
        @model()
        return
    view: () =>
        return [
            Header.getView()
            @table.view()
        ]
