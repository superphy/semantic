class Form
    model: () =>
        @fields ?= []
    constructor: () ->
        @model()
        return
    add_field: (title) =>
        @fields.push new InputGroup(title)
    submit: ->
        if @fields.length > 0
            console.log(input.get_json() for input in @fields)
                
    view: () =>
        return [
            m('form', [
                [input.view(); m('br')] for input in @fields
                m('button[type=button]', {onclick: => @submit()}, 'Submit')
            ])
        ]

class GenomeForm extends Form
    model: () =>
        super()
    constructor:()->
        @model()
        @add_field('Genome Name')
        @add_field('Strain')
        @add_field('Serotype')
        @add_field('Isolation Host')
        @add_field('Isolation Source')
        @add_field('')
    view: ()=>
        #make this one enter parameters for style.
        super()

class InputGroup
    model: (label, initial_value) =>
        @label ?= label
        @value ?= m.prop(initial_value)
    constructor: (
        label='give me a title', initial_value='') ->
        @model(label, initial_value)
    get_json: () =>
        json = {}
        json[@label] = @value()
        return json

    view: () =>
        return [
            m('label', @label)
            m('input', {oninput: m.withAttr('value', @value), value: @value()})
        ]