#pragma once

#include <array>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace almond::fem
{
    struct Node
    {
        double x{0.0};
        double y{0.0};
    };

    struct Element
    {
        std::array<std::size_t, 3> node_ids{};
        double conductivity{1.0};
    };

    class Mesh
    {
    public:
        Mesh() = default;

        Mesh(std::vector<Node> nodes, std::vector<Element> elements)
            : m_nodes(std::move(nodes))
            , m_elements(std::move(elements))
        {
            validate();
        }

        [[nodiscard]] std::span<const Node> nodes() const noexcept { return m_nodes; }
        [[nodiscard]] std::span<Node> nodes() noexcept { return m_nodes; }

        [[nodiscard]] std::span<const Element> elements() const noexcept { return m_elements; }
        [[nodiscard]] std::span<Element> elements() noexcept { return m_elements; }

        void add_node(Node node)
        {
            m_nodes.emplace_back(node);
        }

        void add_element(Element element)
        {
            if (m_nodes.empty())
            {
                throw std::logic_error("Cannot add elements to an empty mesh");
            }

            for (auto id : element.node_ids)
            {
                if (id >= m_nodes.size())
                {
                    throw std::out_of_range("Element references a node that does not exist");
                }
            }

            m_elements.emplace_back(element);
        }

        [[nodiscard]] const Node& node(std::size_t index) const
        {
            if (index >= m_nodes.size())
            {
                throw std::out_of_range("Node index out of range");
            }
            return m_nodes[index];
        }

        [[nodiscard]] const Element& element(std::size_t index) const
        {
            if (index >= m_elements.size())
            {
                throw std::out_of_range("Element index out of range");
            }
            return m_elements[index];
        }

        [[nodiscard]] std::size_t node_count() const noexcept { return m_nodes.size(); }
        [[nodiscard]] std::size_t element_count() const noexcept { return m_elements.size(); }

        void validate() const
        {
            for (const auto& element : m_elements)
            {
                for (auto node_id : element.node_ids)
                {
                    if (node_id >= m_nodes.size())
                    {
                        throw std::out_of_range("Mesh element references an invalid node index");
                    }
                }
            }
        }

    private:
        std::vector<Node> m_nodes{};
        std::vector<Element> m_elements{};
    };
} // namespace almond::fem

