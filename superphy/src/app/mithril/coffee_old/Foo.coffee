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
                    m('button[index=' + index + ']',
                        events.delete_item(@strings),
                        "DELETE"
                    )
                    m('text',string)
                ]) for string, index in @strings
            ])
            #Form to input objects in the array
            m('label', ['New Thing: '
                m('input',
                    { oninput: m.withAttr('value', @value),value : @value()})
                m('button', {
                    disabled : !@value(),
                    onclick: =>
                        if @value() then @strings.push(@value()) and @value("")
                },["Add"])
            ])
        ]