module DelayedEvents
import ..AbstractTypes: EventComponent, ReadoutComponent, RFPulseComponent
import ...Variables: VariableType, get_free_variable, get_pulse, get_readout

struct DelayedEvent{T<:EventComponent} <: EventComponent
    event :: T
    delay :: VariableType
    function DelayedEvent(event::EventComponent, delay=nothing)
        new{typeof(event)}(event, get_free_variable(delay))
    end
end


effective_time(del::DelayedEvent) = effective_time(del.event) + del.delay

get_pulse(del::DelayedEvent{T}) where {T<:RFPulseComponent} = del.event
get_readout(del::DelayedEvent{T}) where {T<:ReadoutComponent} = del.event


end