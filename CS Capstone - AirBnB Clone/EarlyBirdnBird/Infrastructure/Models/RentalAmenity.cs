using System;
using System.Collections.Generic;
using System.ComponentModel.DataAnnotations;
using System.ComponentModel.DataAnnotations.Schema;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


namespace Infrastructure.Models
{
    public class RentalAmenity
    {

        [Key]
        public int Id { get; set; }

        [ForeignKey("RentalId")]
        public int RentalId { get; set; }

        public virtual Rental? Rental { get; set; }

        [ForeignKey("AmenityId")]
        public int AmenityId { get; set; }

        public virtual Amenity? Amenity { get; set; }

    }
}